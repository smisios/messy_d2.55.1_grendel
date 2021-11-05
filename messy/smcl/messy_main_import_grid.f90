! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_MAIN_IMPORT_GRID
! ------------------------------------------------------------------
! ******************************************************************
! Author: Patrick Joeckel, MPICH, Mainz, October 2002
!         extended for curvilinear ocean grids by
!         Astrid  Kerkweg, UniMz, Mainz, 2012-2013
!                          UniBN, Bonn,  2016-2019
!         -> re-structured/expanded for more general GRID application
!                          FZJ-IEK8, Juelich,  2019-
!         -> expanded for pressure Nx2D fields
! ******************************************************************

  USE messy_main_constants_mem, ONLY: DP, STRLEN_MEDIUM, STRLEN_ULONG
  USE messy_main_import_grid_par

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  INTRINSIC :: ABS, ACHAR, ALLOCATED, ASSOCIATED, DATE_AND_TIME, INDEX &
             , LEN_TRIM, MIN, NULL, PRESENT, SIZE, TINY, TRIM, ADJUSTL

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodstr = 'import_grid'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodver = '1.3'

  INTEGER,  PARAMETER            :: NCCNTMAXNLEN =  2*STRLEN_MEDIUM + 1 + 4

  ! TYPE DECLARATION TO HOLD COUNTER INFORMATION
  TYPE t_ncrgcnt
     CHARACTER(LEN=NCCNTMAXNLEN) :: name = ''
     INTEGER                     :: start = 1
     INTEGER                     :: step  = 1
     INTEGER                     :: reset = 1
     INTEGER                     :: current = 1
  END TYPE t_ncrgcnt

  INTEGER,  PUBLIC, SAVE      :: RG_CTRL, RG_NML, RG_STATUS
  INTEGER,  PUBLIC, PARAMETER :: RG_SCAN = 1
  INTEGER,  PUBLIC, PARAMETER :: RG_PROC = 2
  INTEGER,  PUBLIC, PARAMETER :: RG_STOP = 3
  INTEGER,  PUBLIC, PARAMETER :: NML_NEXT = 4
  INTEGER,  PUBLIC, PARAMETER :: NML_STAY = 5
  INTEGER,  PUBLIC, PARAMETER :: RGSTAT_START = 6
  INTEGER,  PUBLIC, PARAMETER :: RGSTAT_CONT = 7
  INTEGER,  PUBLIC, PARAMETER :: RGSTAT_STOP = 8

  PUBLIC :: parse_str

  ! COUNTER INFORMATION HANDLING
  PUBLIC  :: t_ncrgcnt          ! TYPE STRUCT FOR COUNTER INFORMATION I/O
  PUBLIC  :: RGTOOL_NCRGCNT_RST ! MANAGES COUNTER INFORMATION RESTART I/O
  !
  PUBLIC :: PARSE_VARSTR
  PUBLIC :: READ_NAMELIST
  PUBLIC :: READ_CONTROL

  PUBLIC  :: RGTOOL_READ_NCVAR  ! NCREGRID ONE STEP OF ONE FIELD
  PUBLIC  :: RGTOOL_READ_NCFILE ! NCREGRID ONE STEP OF ONE FILE

CONTAINS

  ! ----------------------------------------------------------------------
  SUBROUTINE parse_str(status, strlen, str, nml, var, file, z, heights &
       , p, plevs, ipol, grdname)

    USE messy_main_tools,      ONLY: STRCRACK
    USE messy_main_grid_trafo, ONLY: GTRFCHAR, GTRF_MAX

    IMPLICIT NONE

    INTRINSIC :: TRIM, ADJUSTL

    ! I/O
    INTEGER,          INTENT(OUT)   :: status   ! status information
    INTEGER,          INTENT(IN)    :: strlen   ! max. length of strings
    CHARACTER(LEN=*), INTENT(IN)    :: str      ! string to parse
    CHARACTER(LEN=*), INTENT(INOUT) :: nml      ! namelist file
    CHARACTER(LEN=*), INTENT(INOUT) :: var      ! netCDF variable
    CHARACTER(LEN=*), INTENT(INOUT) :: file     ! netCDF file
    REAL(DP), DIMENSION(:), POINTER :: z        ! level information
    CHARACTER(LEN=*), INTENT(INOUT) :: heights  ! string containing heights
    REAL(DP), DIMENSION(:), POINTER :: p        ! level information
    CHARACTER(LEN=*), INTENT(INOUT) :: plevs    ! string containing heights
    INTEGER         , INTENT(INOUT) :: ipol     ! interpolation method:
                                                ! NRGD : NREGRID
                                                ! SCRP : SCRIP
                                                ! NONE : no regridding
    CHARACTER(LEN=STRLEN_MEDIUM), INTENT(INOUT) :: grdname

    ! LOCAL
    CHARACTER(LEN=strlen), POINTER :: sl1(:)
    CHARACTER(LEN=strlen), POINTER :: sl2(:)
    CHARACTER(LEN=strlen), POINTER :: sl3(:)
    INTEGER                        :: n, m, l
    INTEGER                        :: i, j
    INTEGER                        :: iostat

    status = 1 ! ERROR

    NULLIFY(sl1)
    NULLIFY(sl2)
    NULLIFY(sl3)

    IF (ASSOCIATED(z)) THEN
       DEALLOCATE(z)
       NULLIFY(z)
    END IF

    IF (ASSOCIATED(p)) THEN
       DEALLOCATE(p)
       NULLIFY(p)
    END IF

    CALL STRCRACK(str, ';', sl1, n)
    DO i=1, n

       CALL STRCRACK(sl1(i), '=', sl2, m)
       IF (SIZE(sl2) == 2) THEN
          IF (TRIM(ADJUSTL(sl2(2))) == '') THEN
             status = 2    ! EMPTY SPECIFICATION
             RETURN
          END IF
       END IF

       SELECT CASE(TRIM(ADJUSTL(sl2(1))))
          CASE('NML')
             nml = TRIM(ADJUSTL(sl2(2)))
          CASE('VAR')
             var = TRIM(ADJUSTL(sl2(2)))
          CASE('FILE')
             file = TRIM(ADJUSTL(sl2(2)))
          CASE('Z')
             heights=TRIM(ADJUSTL(sl2(2)))
             CALL STRCRACK(TRIM(ADJUSTL(sl2(2))), ',', sl3, l)
             ALLOCATE(z(l))
             DO j=1, l
                READ(sl3(j),*,IOSTAT=iostat) z(j)
                IF (iostat /= 0) THEN
                   status = 6  ! ERROR IN READING REAL
                   RETURN
                END IF
             END DO
             IF (ASSOCIATED(sl3)) THEN
                DEALLOCATE(sl3) ; NULLIFY(sl3)
             END IF
          CASE('P')
             plevs=TRIM(ADJUSTL(sl2(2)))
             CALL STRCRACK(TRIM(ADJUSTL(sl2(2))), ',', sl3, l)
             ALLOCATE(p(l))
             DO j=1, l
                READ(sl3(j),*,IOSTAT=iostat) p(j)
                IF (iostat /= 0) THEN
                   status = 6  ! ERROR IN READING REAL
                   RETURN
                END IF
             END DO
             IF (ASSOCIATED(sl3)) THEN
                DEALLOCATE(sl3) ; NULLIFY(sl3)
             END IF
          CASE('IPOL')
             DO j=1,GTRF_MAX
                IF (TRIM(ADJUSTL(sl2(2))) == GTRFCHAR(j)) THEN
                   ipol = j
                   EXIT
                ENDIF
             END DO
          CASE ('GRID')
             grdname=TRIM(ADJUSTL(sl2(2)))
          CASE DEFAULT
             write (*,*) 'UNKNOWN SPECIFIER ',sl2(1)
             status = 3 ! UNKNOWN SPECIFIER
             RETURN
       END SELECT

    END DO

    ! CLEAN UP
    IF (ASSOCIATED(sl1)) THEN
       DEALLOCATE(sl1) ; NULLIFY(sl1)
    END IF
    IF (ASSOCIATED(sl2)) THEN
       DEALLOCATE(sl2) ; NULLIFY(sl2)
    END IF

    status = 0 ! NO ERROR

  END SUBROUTINE parse_str
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_NCRGCNT_RST(start, restart, event, c, r, lout, linit)

    ! MANAGES I/O OF COUNTER INFORMATION AT START AND RESTART
    !
    ! Author: Patrick Joeckel, MPICH, September 2005

    USE messy_main_grid_netcdf, ONLY: RGMSG, RGMLI, RGMLIC

    IMPLICIT NONE

    INTRINSIC :: PRESENT, TRIM, REAL

    ! I/O
    LOGICAL,          INTENT(IN)    :: start   ! .true. at first time step
    LOGICAL,          INTENT(IN)    :: restart ! .true. at first time step of
                                               ! rerun
    LOGICAL,          INTENT(IN)    :: event   ! .true. on event
    TYPE(t_ncrgcnt),  INTENT(INOUT) :: c       ! counter - struct
    REAL(DP),         POINTER       :: r       ! INOUT ... current state as real
    LOGICAL,          INTENT(IN)           :: lout
    LOGICAL,          INTENT(IN), OPTIONAL :: linit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_NCRGCNT_RST'
    LOGICAL                     :: zlinit

    ! INIT
    IF (PRESENT(linit)) THEN
       zlinit = linit
    ELSE
       zlinit = .FALSE. ! DEFAULT
    END IF

    ! DO NOTHING IN CASE ...
    IF (.NOT.(start.OR.restart.OR.event.OR.zlinit)) THEN
       IF (lout) CALL RGMSG(substr, RGMLI, &
            'RGTEVENT '''//TRIM(c%name)//''' NOT ACTIVE !')
       RETURN
    ELSE
       IF (lout) CALL RGMSG(substr, RGMLI, &
            'RGTEVENT '''//TRIM(c%name)//''' ACTIVE !')
    END IF

    ! EVENT
    IF (event.AND.(.NOT.(start))) THEN
       IF (lout) CALL RGMSG(substr, RGMLIC, &
            '... UPDATING COUNTER '''//TRIM(c%name)//'''')
       !
       ! UPDATE
       c%current = c%current + c%step                   ! INCREMENT
       IF (c%current > c%reset) c%current = c%start     ! RESET
       !
       IF (ASSOCIATED(r)) r = REAL(c%current,DP)
       !
       IF (lout)  THEN
          CALL RGMSG(substr, RGMLIC, &
               '         START  : ',c%start,' ')
          CALL RGMSG(substr, RGMLIC, &
               '         STEP   : ',c%step,' ')
          CALL RGMSG(substr, RGMLIC, &
               '         RESET  : ',c%reset,' ')
          CALL RGMSG(substr, RGMLIC, &
               '         CURRENT: ',c%current,' ')
       END IF
    END IF

  END SUBROUTINE RGTOOL_NCRGCNT_RST
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  ! ----------------------------------------------------------------------

! ******************************************************************
SUBROUTINE READ_CONTROL(GCTRL, GNML, GSTAT              &
                        ,var                            &
                        ,gridin                         &
                        ,gridout                        &
                        ,RG_TYPE                        &
                        ,lint                           &
                        ,nmlfile                        &
                        ,infostr                        &
                        ,tc, tmin, tmax, tstep, tret    &
                        ,iounit                         &
                        ,ldomainpar                     &
                        ,lvarparallel                   &
                        ,lpresaxis                      &
                        ,lwork                          &
                        ,mnc                            &
                        )

  USE messy_main_constants_mem, ONLY: dp, i8
  USE messy_main_grid_trafo,    ONLY: CHECK_GEOHYBGRID &
                                    , CHECK_NCVAR_ON_GEOHYBGRID &
                                    , REDUCE_INPUT_GRID
  USE messy_main_grid,          ONLY: t_geohybgrid, COPY_GEOHYBGRID           &
                                    , IMPORT_GEOHYBGRID, INIT_GEOHYBGRID
  USE messy_main_grid_netcdf,   ONLY: RGMSG, RGMLE, RGMLEC, RGMLI,  RGMLIC    &
                                    , RGMSG, RGMLW, RGMLWC, RGMLVL, RGMLVM    &
                                    , RGMLVMC, RGMLVLC, RGMSG, ERRMSG         &
                                    , GRD_MAXSTRLEN                           &
                                    , t_ncvar, t_narray, INIT_NARRAY          &
                                    , INIT_NCATT, ADD_NCATT, INIT_NCVAR       &
                                    , IMPORT_NCVAR, RENAME_NCVAR, COPY_NCVAR  &
                                    , SCALE_NARRAY, COPY_NCDIM                &
                                    , NULL_DIMID , nf90_inq_libvers           &
                                    , MAIN_GRID_SET_MESSAGEMODE               &
                                    , NF90_GLOBAL, VTYPE_CHAR                 &
                                    , NF90_CHAR, CAT_NARRAY, MSGMODE_VL       &
                                    , MSGMODE_S, MSGMODE_E                    &
                                    ! KEEP AS COMMENT FOR DEBUGGING
                                    !, MSGMODE_W , MSGMODE_VM , MSGMODE_I     &
                                    , t_ncatt_array, t_multinc, multinetcdf   &
                                    , init_multinc
  USE messy_main_tools,         ONLY: int2str

  IMPLICIT NONE

  ! I/O
  INTEGER, INTENT(IN)                   :: GCTRL   ! RG_SCAN, RG_PROC, RG_STOP
  INTEGER, INTENT(IN)                   :: GNML    ! NML_NEXT, NML_STAY
  INTEGER, INTENT(INOUT)                :: GSTAT   ! RGSTAT_START, RGSTAT_CONT
                                                   ! RGSTAT_STOP
  TYPE (t_ncvar), DIMENSION(:), POINTER :: var     ! list of variables
  TYPE (t_geohybgrid), INTENT(INOUT)    :: gridin  ! input grid info
  TYPE (t_geohybgrid), INTENT(INOUT)    :: gridout ! output grid info
  INTEGER, DIMENSION(:),   POINTER      :: RG_TYPE       ! regridding type
  LOGICAL, INTENT(OUT)                  :: lint    ! input time ?
  CHARACTER(LEN=*) ,    INTENT(IN)      :: nmlfile ! namelist file
  CHARACTER(LEN=STRLEN_ULONG), INTENT(OUT) :: infostr
  !
  ! OPTIONAL I/O
  INTEGER          , INTENT(IN) , OPTIONAL :: tc          ! current step
  INTEGER          , INTENT(OUT), OPTIONAL :: tmin        ! start time step
  INTEGER          , INTENT(OUT), OPTIONAL :: tmax        ! stop time step
  INTEGER          , INTENT(OUT), OPTIONAL :: tstep       ! time step
  INTEGER          , INTENT(OUT), OPTIONAL :: tret        ! return time step
  INTEGER          , INTENT(IN) , OPTIONAL :: iounit      ! I/O unit for nmlfile
  LOGICAL          , INTENT(IN),  OPTIONAL :: ldomainpar  ! parallel application
  LOGICAL          , INTENT(IN),  OPTIONAL :: lvarparallel! parallel application

  ! use pressure axis (only required for NCREGRID)
  LOGICAL          , INTENT(OUT), OPTIONAL :: lpresaxis
  LOGICAL          , INTENT(OUT), OPTIONAL :: lwork

  ! for multi-netcdf descriptor files
  TYPE(t_multinc), INTENT(INOUT), OPTIONAL :: mnc

  ! LOCAL (for DEFAULT)
  CHARACTER(LEN=*), PARAMETER :: substr = 'READ_CONTROL'
  INTEGER                     :: iou   ! I/O unit for namelist file
  ! LOCAL
  LOGICAL                     :: lex   ! file exists ?
  LOGICAL                     :: lopn  ! file opened ?
  INTEGER                     :: iout  ! I/O unit
  INTEGER                     :: nvars ! number of variables
  INTEGER                     :: mvars ! 'effective' number of variables
  INTEGER                     :: status
  TYPE (t_ncvar), DIMENSION(:), POINTER  :: xivar ! list of variables (input)
  TYPE (t_ncvar)                         :: tvar  ! temporal variable

  ! NAMELIST
  INTEGER, SAVE :: fstat ! status of namelist file
  TYPE (t_geohybgrid)                     :: gi      ! input grid
  TYPE (t_geohybgrid),               SAVE :: gi0     ! input grid (raw)
  TYPE (t_geohybgrid)                     :: gg      ! grdfile grid
  TYPE (t_geohybgrid),               SAVE :: gg0     ! grdfile grid
  CHARACTER (LEN=100*GRD_MAXSTRLEN), SAVE :: varstr  ! variable string
  CHARACTER (LEN=GRD_MAXSTRLEN),     SAVE :: outfile ! netCDF output-file
  INTEGER,                           SAVE :: i_t(4)  ! input time control
  INTEGER,                           SAVE :: g_t(3)  ! grid time control
  INTEGER,                           SAVE :: o_t(3)  ! output time control
  LOGICAL,                           SAVE :: lp      ! pressure or sigma
  TYPE (t_narray),                   SAVE :: nmlstr  ! namelist as string
  !
  ! TIME CONTROL
  INTEGER, SAVE                           :: ttc
  INTEGER, SAVE                           :: ttmin, ttmax, ttstep, ttret
  INTEGER                                 :: tg       ! grid time step
  INTEGER                                 :: to       ! output time step
  LOGICAL                                 :: ok       ! result OK ?
  LOGICAL                                 :: ldompar = .FALSE.
  LOGICAL                                 :: lvarpar = .FALSE.
  ! FORCE UNLIMITED ID
  INTEGER                                 :: setuid

  ! VARIABLES
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER, SAVE :: ivar(:) ! input variable names
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER   :: ovar(:)   ! output variable names
  REAL(dp)                    , POINTER   :: scl(:)    ! scaling
  ! range for ixf regridding
  INTEGER                     , POINTER   :: rg_range(:,:)
  ! namelist assign attributes
  TYPE(t_ncatt_array), DIMENSION(:), POINTER  :: namatts => NULL()

  INTEGER                     , POINTER   :: RGT(:)    ! regridding type
  CHARACTER (LEN=3)           , POINTER   :: RGTstr(:) ! regridding type string
  ! order of grid dimensions
  INTEGER, DIMENSION(:)       , POINTER   :: dims
  INTEGER                                 :: axes(3)   ! grid axes of variable
  INTEGER, DIMENSION(2)                   :: zstart, zcount
  CHARACTER(LEN=GRD_MAXSTRLEN)            :: att_c= ' '

  ! multi-netcdf file descriptor option currently only implemented for infile;
  ! (gridfile could be trated similarly, outfile currently not possible,
  !  since number of time steps are read from netcdf files, see
  !  subroutine multinetcdf)
  TYPE(t_multinc)              :: mcin     ! multi-netcdf descriptor
  CHARACTER(LEN=GRD_MAXSTRLEN) :: fnamenc  ! netcdf filename
  INTEGER                      :: ztf      ! time step number in file

  IF (PRESENT(ldomainpar)) THEN
     ldompar =ldomainpar
  ELSE
     ldompar = .FALSE.
  ENDIF
  IF (PRESENT(lvarparallel)) THEN
     lvarpar = lvarparallel
  ELSE
     lvarpar = .FALSE.
  ENDIF
  ! parallelise only once (domain  => ldompar = .TRUE.
  !                     or species => lvarpar = .TRUE. )
  IF (ldompar) lvarpar = .FALSE.

  IF(lvarpar .OR. ldompar)  THEN
     me    = my_mpi%rank
     nproc = my_mpi%nproc
     comm  = my_mpi%comm
  ELSE
     me    = 0
     nproc = 1
     comm  = -1
  END IF

  CALL READ_CONTROL_INIT

  IF (i_am_worker) THEN
     CALL READ_CONTROL_WORK
  END IF

  IF (PRESENT(lwork)) lwork = i_am_worker

  RETURN

 CONTAINS

 SUBROUTINE READ_CONTROL_INIT

   USE MESSY_MAIN_TOOLS, ONLY: position_nml &
                             , POSITIONED, MISSING, LENGTH_ERROR, READ_ERROR

   IMPLICIT NONE

   INTEGER :: i, j
   INTEGER, SAVE :: n
   INTEGER :: ierr

  ! INIT
  CALL MAIN_GRID_SET_MESSAGEMODE(MSGMODE_S + MSGMODE_E + MSGMODE_VL) ! &!)
                                ! + MSGMODE_W + MSGMODE_VM + MSGMODE_I)
  ! CHECK INTERFACE
  IF ((GNML /= NML_NEXT).AND.(GNML /= NML_STAY)) THEN
     CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED VALUE OF RG_NML !')
  END IF

  ! SET DEFAULTS
  IF (PRESENT(iounit)) THEN
     iou = iounit
  ELSE
     iou = 17                  ! DEFAULT
  END IF
  !
  ! INITIALIZE OUTPUT PARAMETER
  IF (ASSOCIATED(var)) THEN
     DO i=1, SIZE(var)
        CALL INIT_NCVAR(var(i))
     END DO
     DEALLOCATE(var, STAT=status)
     CALL ERRMSG(substr,status,1)
  END IF
  NULLIFY(var)
  CALL INIT_GEOHYBGRID(gridin)

  ! INIT
  NULLIFY(dims)
  NULLIFY(ivar)
  NULLIFY(ovar)
  NULLIFY(scl)
  NULLIFY(rg_range)

  IF (ASSOCIATED(namatts)) THEN
     DO i = 1, SIZE(namatts)
        DO j = 1, namatts(i)%natts
           CALL INIT_NCATT(namatts(i)%atts(j))
        END DO
        NULLIFY(namatts(i)%atts)
     END DO
     DEALLOCATE(namatts)
  END IF
  NULLIFY(namatts)

  NULLIFY(RGT)
  NULLIFY(RGTstr)
  NULLIFY(xivar)

  ! CHECK NAMELIST FILE; OPEN FILE
  INQUIRE(file=TRIM(nmlfile), exist=lex, opened=lopn, number=iout)
  IF (.NOT.lex) THEN
     CALL RGMSG(substr, RGMLE, &
          'NAMELIST-FILE '''//TRIM(nmlfile)//''' NOT FOUND !')
  END IF
  IF (lopn) THEN  ! FILE OPEN
     IF (iout /= iou) THEN   ! CHECK UNIT
        CALL RGMSG(substr, RGMLE, &
             'NAMELIST-FILE '''//TRIM(nmlfile)//''' ALREADY OPEN ON UNIT', &
             iout,' !')
     ELSE
        GSTAT = RGSTAT_CONT
     END IF
  ELSE            ! FILE NOT OPEN
     !
     n = 0 ! intialise for n'th namelist read
     !
     CALL RGMSG(substr, RGMLVL, &
         '-------------------------------------------------------')
     CALL RGMSG(substr, RGMLVL, 'START READ PROCEDURE FOR REGRIDDING ')
     CALL RGMSG(substr, RGMLVL, &
          ' <Author: Patrick Joeckel, MPICH, June 2002>' )
     CALL RGMSG(substr, RGMLVL, &
          ' ( COMPILED WITH netCDF Fortran90 LIBRARY' )
     CALL RGMSG(substr, RGMLVL, &
          '   VERSION '//TRIM(NF90_inq_libvers())//' )' )
     CALL RGMSG(substr, RGMLVL, &
          '-------------------------------------------------------')
     CALL RGMSG(substr, RGMLVL, &
          'OPENING NAMELIST FILE '''//TRIM(nmlfile)//''' ON UNIT', &
          iou,' !')
     OPEN(iou,file=TRIM(nmlfile))
     GSTAT = RGSTAT_START
  END IF

  ! READ NEXT NAMELIST
  IF ((GSTAT == RGSTAT_START).OR.(GNML == NML_NEXT)) THEN

     CALL INIT_GEOHYBGRID(gi0)
     CALL INIT_GEOHYBGRID(gg0)
     CALL INIT_NARRAY(nmlstr)

     CALL RGMSG(substr, RGMLI, &
          'READING NAMELIST FROM FILE '''//TRIM(nmlfile)//''' ...')

     CALL position_nml ('regrid', iou, ierr, n)
     SELECT CASE(ierr)
     CASE(POSITIONED)
        ! NOTHING TO DO: continue below to read namelist
     CASE(MISSING)
        !
        i_am_worker = .FALSE.
        IF (n == 1) THEN
           CALL RGMSG(substr, RGMLE, &
                'NO REGRID NAMELIST FOUND IN FILE '''//&
                &TRIM(nmlfile)//''' !' , .TRUE.)
        ELSE
           CALL RGMSG(substr, RGMLVL, &
                'END OF FILE '''//TRIM(nmlfile)//''' REACHED !' )
        ENDIF
        GSTAT = RGSTAT_STOP
        RETURN ! qqq?
     CASE(READ_ERROR)
        CALL RGMSG(substr, RGMLE, &
             'READ ERROR READING NAMELIST FROM FILE '''//&
             &TRIM(nmlfile)//''' !' , .TRUE.)
     CASE(LENGTH_ERROR)
        CALL RGMSG(substr, RGMLE, &
                'LENGTH ERROR READING NAMELIST FROM FILE '''//&
                &TRIM(nmlfile)//''' !' , .TRUE.)
     END SELECT

     CALL READ_NAMELIST(iou, fstat               &
          ,gi0, gg0, outfile       &
          ,varstr                  &
          ,i_t, g_t, o_t, lp, lint, nmlstr &
          ,infostr)

     IF (fstat == 0) THEN
        ! --- TIME CONTROL ---------------------------------
        ! UPDATE TIME CONTROL FROM NAMELIST
        ttmin = i_t(1)
        ttstep = i_t(2)
        ttmax = i_t(3)
        ttret = i_t(4)
        !
        ! OVERWRITE TIME CONTROL FROM SUBROUTINE CALL
        IF (PRESENT(tmin)) THEN
           tmin = ttmin
        END IF
        IF (PRESENT(tmax)) THEN
           tmax = ttmax
        END IF
        IF (PRESENT(tstep)) THEN
           tstep = ttstep
        END IF
        IF (PRESENT(tret)) THEN
           tret = ttret
        END IF
        !
        ! CHECK TIME CONTROL SETTINGS
        ! CHECK INFILE TIME STEPPING / INTERFACE MODE TIME STEPPING
        IF (ttmin <= 0) THEN
           CALL RGMSG(substr, RGMLE,  'i_t(1) IN NAMELIST OR tmin AT ', .false.)
           CALL RGMSG(substr, RGMLEC, 'SUBROUTINE CALL MUST BE >0 !')
        END IF
        IF (ttstep <= 0) THEN
           CALL RGMSG(substr, RGMLE,  'i_t(2) IN NAMELIST OR tstep AT ',.false.)
           CALL RGMSG(substr, RGMLEC, 'SUBROUTINE CALL MUST BE >0 !')
        END IF
        IF (ttmax < 0) THEN
           CALL RGMSG(substr, RGMLE,  'i_t(3) IN NAMELIST OR tmax AT ', .false.)
           CALL RGMSG(substr, RGMLEC, 'SUBROUTINE CALL MUST BE >=0 !')
        END IF
        IF (ttret < 0) THEN
           CALL RGMSG(substr, RGMLE,  'i_t(4) IN NAMELIST OR tret AT ', .false.)
           CALL RGMSG(substr, RGMLEC, 'SUBROUTINE CALL MUST BE >=0 !')
        END IF
        !
        IF (ttmax == 0) THEN  ! SPECIAL CASE
           CALL RGMSG(substr, RGMLW,  'tmax (i_t(3)) RESET TO tmin (i_t(1)) !')
           ttmax = ttmin
        END IF
        IF (ttret == 0) THEN  ! SPECIAL CASE
           CALL RGMSG(substr, RGMLW,  'tret (i_t(4)) RESET TO tmin (i_t(1)) !')
           ttret = ttmin
        END IF
        !
        IF (ttmax < ttmin) THEN
           CALL RGMSG(substr, RGMLE,  'i_t(3) IN NAMELIST OR tmax AT ', .false.)
           CALL RGMSG(substr, RGMLEC &
                , 'SUBROUTINE CALL MUST BE >= i_t(1) OR tmin!')
        END IF
        IF ( .NOT.( (ttret >= ttmin).AND.(ttret <= ttmax) ) ) THEN
           CALL RGMSG(substr, RGMLE, &
                'NCREGRID IN INTERFACE MODE ONLY RETURNS DATA, IF', .false.)
           CALL RGMSG(substr, RGMLEC, &
                '   i_t(1) <= i_t(4) <= i_t(3)  IN THE NAMELIST, OR', .false.)
           CALL RGMSG(substr, RGMLEC, &
                '   tmin   <= tret   <= tmax    AT SUBROUTINE CALL!', .false.)
           CALL RGMSG(substr, RGMLEC, &
                'PLEASE CORRECT NAMELIST (MOST PROBABLY i_t(4)), OR', .false.)
           CALL RGMSG(substr, RGMLEC, &
                'PARAMETERS IN SUBROUTINE CALL (tret) ACCORDINGLY!')
        END IF
        ! CHECK GRDFILE TIME STEPPING
        ! NOTE: IN PRINCIPLE IT IS POSSIBLE TO SPECIFY A GRDFILE
        !       (IN THE NAMELIST)
        !       IN INTERFACE MODE; THIS OVERWRITES THE GRID SPECIFIED BY THE
        !       INTERFACE.
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
        ! CHECK OUTFILE TIME STEPPING
        IF (o_t(1) <= 0) THEN
           CALL RGMSG(substr, RGMLE, 'o_t(1) IN NAMELIST MUST BE >0 !')
        END IF
        IF (o_t(1) <= 0) THEN
           CALL RGMSG(substr, RGMLE, 'o_t(2) IN NAMELIST MUST BE >0 !')
        END IF
        ! o_t(3) IS NOT USED
        ! --- END TIME CONTROL ---------------------------------

        IF (PRESENT(tc)) THEN
           ttc = tc
        ELSE
           ttc = ttmin
        ENDIF

     END IF

     ! CHECK END OF FILE
     IF (fstat == 0) THEN    ! GO ON
        CALL RGMSG(substr, RGMLIC, '... OK !')

        ! CHECK MULTI-NETCDF DESCRIPTOR
        IF (PRESENT(mnc)) THEN
           CALL multinetcdf(status, gi0%file, gi0%timem%name &
                , fnamenc, ttc, ztf, mnc)
        ELSE
           CALL multinetcdf(status, gi0%file, gi0%timem%name &
                , fnamenc, ttc, ztf, mcin)
        END IF
        IF (status /= 0) &
             CALL RGMSG(substr, RGMLE, '... error in multinetcdf !')
        gi0%file = TRIM(fnamenc)
        gi0%t    = ztf
        ttc      = ztf
        IF (.NOT. PRESENT(mnc)) &
             CALL init_multinc(mcin) ! free memory, if local instance was used

        ! GET VARIABLE INFORMATION
        CALL PARSE_VARSTR(varstr, ivar, ovar, scl, namatts   &
             ,RGT, RGTstr, rg_range, TRIM(gi0%file))

        IF (lvarpar .OR. ldompar) THEN
        ! NOTE: THE FOLLOWING SUBROUTINE CALLS ARE REQUIRED TO SET
        !       i_am_worker CORRECTLY
           IF (ldompar) THEN
              i_am_worker = .TRUE.
              ! in parallel domain decomposition mode, each PE must work on
              ! all variables; nproc_work is the number of PEs among which the
              ! number of variables is distributet; thus nproc_work = 1 here.
              nproc_work  = 1
           ELSE
              ! nproc_work and i_am_worker ist set in this subroutine
              CALL DISTRIBUTE_VARS_ON_PES(ivar, ovar, scl &
                   , namatts, RGT, RGTstr, rg_range)
              IF (i_am_worker) THEN
                 DO i=1, SIZE(ivar)
                    CALL RGMSG(substr, RGMLIC &
                         , '... ... '//TRIM(ivar(i))//' on PE=', me, '')
                 END DO
              END IF
           ENDIF
        ELSE
           ! ... in this case, the routine has been called only on p_io ...
           i_am_worker = .TRUE.
           nproc_work  = 1
        ENDIF
        !
        GSTAT = RGSTAT_CONT
        IF (PRESENT(lpresaxis)) lpresaxis = lp
        !
     ELSE                    ! ERROR OR END OF FILE
        i_am_worker = .FALSE.
        IF (fstat > 0) THEN  ! ERROR
           CLOSE(iou)
           CALL RGMSG(substr, RGMLE, &
                'READ ERROR READING NAMELIST FROM FILE '''//&
                &TRIM(nmlfile)//''' !' , .false.)
           CALL RGMSG(substr, RGMLEC, &
                'STATUS: ',fstat,' ')
        ELSE                 ! END OF FILE
           CALL RGMSG(substr, RGMLVL, &
                'END OF FILE '''//TRIM(nmlfile)//''' REACHED !' )
        END IF               ! ERROR OR END OF FILE
        GSTAT = RGSTAT_STOP
     END IF
  ELSE IF (GNML == NML_STAY .AND. GCTRL /= RG_STOP) THEN

     ! GET VARIABLE INFORMATION
     CALL PARSE_VARSTR(varstr, ivar, ovar, scl, namatts     &
          ,RGT, RGTstr, rg_range, TRIM(gi0%file))
  END IF

  IF (GCTRL == RG_STOP) GSTAT = RGSTAT_STOP

  ! --------------------------------------------------
  IF (GSTAT == RGSTAT_STOP) THEN

     CLOSE(iou)
     ! CLEAN
     CALL RGMSG(substr, RGMLVM, 'CLEANING MEMORY ...')
     CALL INIT_GEOHYBGRID(gridin)
     !
     IF (ASSOCIATED(var)) THEN
        DO i=1, SIZE(var)
           CALL INIT_NCVAR(var(i))
        END DO
        DEALLOCATE(var, STAT=status)
        CALL ERRMSG(substr,status,2)
     END IF
     IF (ASSOCIATED(RG_TYPE)) THEN
        DEALLOCATE(RG_TYPE)
        NULLIFY(RG_TYPE)
     ENDIF
     NULLIFY(var)
     !
     CALL INIT_GEOHYBGRID(gi0)
     CALL INIT_GEOHYBGRID(gg0)
     CALL INIT_NARRAY(nmlstr)
     !
     tg = 0
     to = 0
     CALL RGMSG(substr, RGMLVM, '... DONE (CLEANING MEMORY) !')
     !
     ! MESSAGE
     CALL RGMSG(substr, RGMLVL, &
         '-------------------------------------------------------')
     CALL RGMSG(substr, RGMLVL, &
          'END REGRIDDING PROCEDURE')
     CALL RGMSG(substr, RGMLVL, &
         '-------------------------------------------------------')
     RETURN
  END IF
  ! --------------------------------------------------
  RETURN

END SUBROUTINE READ_CONTROL_INIT

SUBROUTINE READ_CONTROL_WORK

  IMPLICIT NONE

  INTEGER          :: i, ix, j, jx
  CHARACTER(LEN=5) :: ixfstr(2)

  IF (GSTAT == RGSTAT_STOP) RETURN

  IF (GNML == NML_NEXT .AND. (ldompar .OR. lvarpar)) THEN
     CALL EXPAND_FILENAME(outfile, me)
  END IF

  CALL RGMSG(substr, RGMLIC, '... OUTPUT FILE '//TRIM(outfile))
  CALL RGMSG(substr, RGMLIC, '... OK !')

  nvars = SIZE(ivar)
  ALLOCATE(xivar(nvars), STAT=status)
  CALL ERRMSG(substr,status,3)

!!$  ! START TIME LOOP
!!$  tg = g_t(1) - g_t(2)              ! INITIALIZE time step in grid file
!!$  to = o_t(1) - o_t(2)              ! INITIALIZE time step in output file

  CALL RGMSG(substr, RGMLVL, &
       '-------------------------------------------------------')

  ! copy raw grid from namelist input
  CALL COPY_GEOHYBGRID(gi, gi0)
  CALL COPY_GEOHYBGRID(gg, gg0)

  ! UPDATE OUTPUT TIME STEP
  tg = tg + g_t(2)              ! NEXT TIME STEP IN GRID FILE
  to = to + o_t(2)              ! NEXT OUTPUT TIME STEP
  IF ((tg > g_t(3)).AND.(g_t(3) > 0)) THEN  ! RESET
     tg = g_t(1)
  END IF

  gg%t = tg
  CALL RGMSG(substr, RGMLVM, 'INFILE  TIME STEP: ',gi%t,' ')
  CALL RGMSG(substr, RGMLVM, 'GRDFILE TIME STEP: ',gg%t,' ')
  CALL RGMSG(substr, RGMLVM, 'OUTFILE TIME STEP: ',to,' ')
  IF (lint) THEN
     CALL RGMSG(substr, RGMLVM, ' USING INFILE TIME FOR OUTPUT !')
  ELSE
     CALL RGMSG(substr, RGMLVM, ' USING GRDFILE TIME FOR OUTPUT !')
  END IF

  ! GET INPUT GRID
  CALL RGMSG(substr, RGMLVL, &
       'IMPORTING INPUT GRID ... ')
  CALL RGMSG(substr, RGMLVL, &
       '... IMPORTING FROM '''//TRIM(gi%file)//'''')
  CALL IMPORT_GEOHYBGRID(gi)
  !
  CALL RGMSG(substr, RGMLIC, ' ... CHECKING ...')
  CALL CHECK_GEOHYBGRID(gi)
  CALL RGMSG(substr, RGMLIC, ' ... O.K.')

  ! GET DESTINATION GRID
  CALL RGMSG(substr, RGMLVL, 'IMPORTING OUTPUT GRID ...')
  IF (gridout%ID == -99) THEN
     CALL RGMSG(substr, RGMLVL, &
          '... IMPORTING FROM '''//TRIM(gg%file)//'''')
     CALL IMPORT_GEOHYBGRID(gg)
     CALL COPY_GEOHYBGRID(gridout, gg)
  END IF

  IF (lint) THEN  ! ADJUST OUTPUT TIME STEPPING FOR INTERFACE MODE
     to = ttc
     gridout%t = gi%t
  ELSE
     to = tg
  END IF

  ! REDUCE INPUT GRID ACCORDING TO OUTPUT GRID
  CALL RGMSG(substr, RGMLVL, 'REDUCING INPUT GRID ...')
  CALL REDUCE_INPUT_GRID(gi, gridout)

  ! GET INPUT VAR
  CALL RGMSG(substr, RGMLVL, &
       'IMPORTING ',nvars,' VARIABLE(S) FROM FILE '''//TRIM(gi%file)//'''')
  CALL RGMSG(substr, RGMLVL, &
       '('''//TRIM(gi%timem%name)//''' STEP: ',ttc,')')
  mvars = 0

  DO i=1, nvars  ! LOOP OVER VARIABLE NAMES
     CALL RGMSG(substr, RGMLVL,' ... '''//TRIM(ivar(i))//''' ...')
     IF (ASSOCIATED(gi%timem%dim)) THEN
        setuid = gi%timem%dim(1)%id
     ELSE
        setuid = NULL_DIMID
     END IF
     IF (gi%start(1) > 0 .OR. gi%start(2) > 0 .OR. &
          gi%count(1) > 0 .OR. gi%count(2) > 0)  THEN
        zstart = gi%start
        zcount = gi%count
        CALL IMPORT_NCVAR(tvar, ttc, varname=TRIM(ivar(i)) &
             ,file=TRIM(gi%file)             &
             ,setuid=setuid                  &
             ,pstart=zstart, pcount=zcount, hdimids=gi%hdimids)
     ELSE
        CALL IMPORT_NCVAR(tvar, ttc, varname=TRIM(ivar(i)) &
             ,file=TRIM(gi%file)             &
             ,setuid=setuid                  )
     ENDIF

     ! CHECK
     CALL RGMSG(substr, RGMLIC, '  ... CHECKING FOR GRID CONFORMITY ...')
     CALL CHECK_NCVAR_ON_GEOHYBGRID(tvar, gi, dims, axes, ok)

     ! NOT NEEDED HERE
     axes(:) = 0
     DEALLOCATE(dims, STAT=status)
     CALL ERRMSG(substr,status,4)
     NULLIFY(dims)
     !
     IF (ok) THEN
        ! FILENAME
        CALL ADD_NCATT(tvar, 'mmig_input_file', vs=TRIM(gi%file))
        ! RENAME
        IF (TRIM(ivar(i)) /= TRIM(ovar(i))) THEN
           CALL RGMSG(substr, RGMLVLC, &
                '  ... RENAMING TO '''//TRIM(ovar(i))//'''... ')
           CALL RENAME_NCVAR(tvar, TRIM(ovar(i)))
        END IF
        CALL ADD_NCATT(tvar, 'mmig_variable_name', vs=TRIM(ivar(i)))
        ! SCALE
        IF (ABS(scl(i) - 1.0_dp) >= TINY(1.0_dp)) THEN
           CALL  RGMSG(substr, RGMLVLC,'  ... SCALING BY ',scl(i),' ...')
           CALL SCALE_NARRAY(tvar%dat, scl(i))
           CALL ADD_NCATT(tvar, 'mmig_scaled_by',      &
                vd=(/ scl(i) /))
        END IF
        IF (ANY(rg_range(i,:) /= 0)) THEN
           CALL INT2STR(ixfstr(1), rg_range(i,1))
           CALL INT2STR(ixfstr(2), rg_range(i,2))
           CALL RGMSG(substr, RGMLVLC &
                ,'  ... SETTING IXF RANGE ('//TRIM(ixfstr(1))//&
                &':'//TRIM(ixfstr(2))//') ...')
           CALL ADD_NCATT(tvar, 'mmig_ixf_range' &
                , vi=(/ INT(rg_range(i,1),i8), INT(rg_range(i,2),i8) /))
        END IF
        DO ix = 1, namatts(i)%natts
           att_c =' '
           DO jx = 1,namatts(i)%atts(ix)%dat%n
              att_c(jx:jx) = namatts(i)%atts(ix)%dat%vc(jx)
           END DO
           CALL  RGMSG(substr, RGMLVLC, &
              '  ... assign attribute '//TRIM(namatts(i)%atts(ix)%name)//' ...')
           CALL  RGMSG(substr, RGMLVLC, &
                ' ... with value '//TRIM(att_c)//' ...')
           CALL ADD_NCATT(tvar,TRIM(namatts(i)%atts(ix)%name) &
                   , vs=TRIM(att_c), replace=.TRUE. )
        END DO

        ! ADJUST TIME AXIS
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
              CALL RGMSG(substr, RGMLVMC,'      '//TRIM(tvar%dim(j)%name)&
                   //' (ID = ',j,') -> ')
              CALL RGMSG(substr, RGMLVMC,'      '&
                   &//TRIM(gridout%timem%dim(1)%name)&
                   //' (ID = ',setuid,') ...')
              CALL COPY_NCDIM(tvar%dim(j), gridout%timem%dim(1))
           END IF
        END IF
        ! SAVE IN LIST
        mvars = mvars + 1
        CALL COPY_NCVAR(xivar(mvars), tvar)

        RGT(mvars)    = RGT(i)      ! shift
        RGTstr(mvars) = RGTstr(i)   ! shift
     ELSE
        CALL  RGMSG(substr, RGMLVMC, &
             '  ... NOT CONSISTENT WITH GRID! IGNORING !')
     END IF
     ! CLEAN
     CALL INIT_NCVAR(tvar)
  END DO        ! LOOP OVER VARIABLE NAMES

  ! WHAT TO DO WITH IMPORTED VARIABLES
  ! ***************************************************
  SELECT CASE(GCTRL)
     ! ***************************************************
     ! ---------------------------------------------------
  CASE(RG_SCAN)                ! RETURN INPUT VARIABLES
     ! ---------------------------------------------------
     CALL RGMSG(substr, RGMLVL, 'SCAN MODE ...')
     ! RETURN VAR
     CALL RGMSG(substr, RGMLVM, 'TIME STEP ...',ttret,' ')
     CALL RGMSG(substr, RGMLVM, '... RETURNING RAW DATA')
     ALLOCATE(var(mvars), STAT=status)
     CALL ERRMSG(substr,status,5)
     ALLOCATE(RG_TYPE(mvars), STAT=status)
     CALL ERRMSG(substr,status,6)

     DO i=1, mvars
        CALL INIT_NCVAR(var(i))
        CALL COPY_NCVAR(var(i),xivar(i))
        RG_TYPE(i) = RGT(i)
     END DO
     ! RETURN INPUT GRID
     CALL RGMSG(substr, RGMLVM, '... RETURNING INPUT GRID')
     CALL COPY_GEOHYBGRID(gridin, gi)
     CALL RGMSG(substr, RGMLVL, '... DONE (SCAN MODE) !')

     CALL RGMSG(substr, RGMLVL, 'SETTING OUTFILE DATA (SCAN MODE) !')
     IF (TRIM(outfile) /= '') THEN
        gridout%file = TRIM(outfile)    ! name of output-file
        gridout%t = to                  ! output time step

        ! GET OLD ATTRIBUTE
        CALL CAT_NARRAY(gridout%att%dat, nmlstr)
        gridout%att%dat%type = VTYPE_CHAR
        gridout%att%xtype   = NF90_CHAR
        gridout%att%len     = gridout%att%dat%dim(1)
        gridout%att%name    = 'RG_HISTORY'
        gridout%att%varid   = NF90_GLOBAL
     END IF

     CALL RGMSG(substr, RGMLVL, ' DONE SETTING OUTFILE DATA (SCAN MODE) !')

     ! ---------------------------------------------------
  CASE(RG_PROC)                ! REGRID
     ! NOTHING TO DO in READ_CONTROL
     ! ---------------------------------------------------
     ! ---------------------------------------------------
  CASE (RG_STOP)
     ! DO NOTHING; THIS POINT IS NOT REACHED
     ! ---------------------------------------------------
     ! ---------------------------------------------------
  CASE DEFAULT
     ! ---------------------------------------------------
     CALL RGMSG(substr, RGMLE, 'VALUE OF RG_CTRL NOT RECOGNIZED !')
     ! ***************************************************
  END SELECT
     ! ***************************************************

  ! CLEAN
  IF (ASSOCIATED(xivar)) THEN
     DO i=1, SIZE(xivar)
        CALL INIT_NCVAR(xivar(i))
     END DO
     ! DO NOT DEALLOCATE/NULLIFY 'xivar' HERE, SINCE IT
     ! HAS BEEN ALLOCATED OUTSIDE THE TIME LOOP !
  END IF
  !
  CALL INIT_GEOHYBGRID(gg)
  CALL INIT_GEOHYBGRID(gi)

  ! CLEAN
  IF (ASSOCIATED(ivar)) DEALLOCATE(ivar)
  IF (ASSOCIATED(ovar)) DEALLOCATE(ovar)
  IF (ASSOCIATED(scl))  DEALLOCATE(scl)
  IF (ASSOCIATED(rg_range))  DEALLOCATE(rg_range)

  IF (ASSOCIATED(namatts)) THEN
     DO i = 1, SIZE(namatts)
        DO j = 1, namatts(i)%natts
           CALL INIT_NCATT(namatts(i)%atts(j))
        END DO
        NULLIFY(namatts(i)%atts)
     END DO
     DEALLOCATE(namatts)
  END IF
  NULLIFY(namatts)

  IF (ASSOCIATED(RGT))  DEALLOCATE(RGT)
  IF (ASSOCIATED(RGTstr)) DEALLOCATE(RGTstr)
  NULLIFY(ivar, ovar, scl, RGT, RGTstr, rg_range)

  IF (ASSOCIATED(xivar)) THEN
     ! CALL INIT_NCVAR(xivar(i)) HAS ALREADY BEEN DONE INSIDE
     ! (AT THE END OF) THE TIME LOOP
     DEALLOCATE(xivar, STAT=status)
     CALL ERRMSG(substr,status,7)
     NULLIFY(xivar)
  END IF
  CALL RGMSG(substr, RGMLVL, ' DONE  READ_CONTROL_WORK !')

  END SUBROUTINE READ_CONTROL_WORK

END SUBROUTINE READ_CONTROL


! ******************************************************************
! ------------------------------------------------------------------
SUBROUTINE PARSE_VARSTR(var, ivar, ovar, scl, namatts, RGT, RGTstr &
     , RGrange, file)

  USE messy_main_grid_netcdf,  ONLY: GRD_MAXSTRLEN, t_ncvar, ERRMSG            &
                                   , RGMSG, RGMLI, RGMLE, RGMLVM, RGMLW, RGMLWC&
                                   , SCAN_NCVAR, INIT_NCVAR, IMPORT_NCVAR      &
                                   , STRING, INIT_NARRAY, t_ncatt_array        &
                                   , VTYPE_CHAR, INIT_NCATT
  USE messy_main_grid_trafo,   ONLY: RG_INT, RG_EXT, RG_IDX, RG_IXF

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)           :: var       ! variable string
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER  :: ivar(:)   ! input variable names
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER  :: ovar(:)   ! output variabel names
  REAL(dp)                    , POINTER  :: scl(:)    ! scaling factor
  TYPE (t_ncatt_array)        , POINTER  :: namatts(:) ! array of attributes
  INTEGER                     , POINTER  :: RGT(:)    ! regridding type
  CHARACTER(LEN=3)            , POINTER  :: RGTstr(:) ! regridding type string
  CHARACTER(LEN=21)                      :: RGTrangestr
  ! regridding range (for IXF)
  INTEGER                     , POINTER  :: RGrange(:,:)
  CHARACTER(LEN=*)            , OPTIONAL, INTENT(IN) :: file ! filename

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER           :: substr = 'PARSE_VARSTR'
  TYPE (t_ncvar), DIMENSION(:), POINTER :: kvar   ! variables from file
  INTEGER                               :: nvar   ! number of variables
  INTEGER                               :: status
  INTEGER                               :: i,j,ix, jx
  INTEGER                               :: idx1, idx2, idx3, idx4
  INTEGER                               :: idxl, idxr
  CHARACTER(GRD_MAXSTRLEN), ALLOCATABLE :: lvar(:)   ! LIST OF VARIABLE INFOs
  ! LIST OF VAR attribute INFOs
  INTEGER                               :: natt
  INTEGER                               :: idxatt
  CHARACTER(GRD_MAXSTRLEN), ALLOCATABLE :: attvar(:)
  CHARACTER(GRD_MAXSTRLEN), ALLOCATABLE :: varatt(:)
  CHARACTER(GRD_MAXSTRLEN)              :: str = ' '
  INTEGER                               :: len
  CHARACTER(GRD_MAXSTRLEN)              :: infile    ! filename

  ! INIT
  IF (PRESENT(file)) THEN
     infile = TRIM(file)
  ELSE
     infile = ''
  END IF
  NULLIFY(kvar)

  IF (ASSOCIATED(ivar))   DEALLOCATE(ivar);   NULLIFY(ivar)
  IF (ASSOCIATED(ovar))   DEALLOCATE(ovar);   NULLIFY(ovar)
  IF (ASSOCIATED(scl))    DEALLOCATE(scl);    NULLIFY(scl)
  IF (ASSOCIATED(RGT))    DEALLOCATE(RGT);    NULLIFY(RGT)
  IF (ASSOCIATED(RGTstr)) DEALLOCATE(RGTstr); NULLIFY(RGTstr)
  IF (ASSOCIATED(RGrange)) DEALLOCATE(RGrange); NULLIFY(RGrange)

  IF (ASSOCIATED(namatts)) THEN
     DO i = 1, SIZE(namatts)
        DO j = 1, namatts(i)%natts
           CALL INIT_NCATT(namatts(i)%atts(j))
        END DO
        DEALLOCATE(namatts(i)%atts); NULLIFY(namatts(i)%atts)
     END DO
     DEALLOCATE(namatts)
  END IF
  NULLIFY(namatts)

  !  count number of variables
  IF ((TRIM(var) == '').OR.(TRIM(var) == ';')) THEN   ! EMPTY VAR LIST
     ! NO VARIABLES GIVEN -> SCAN FOR ALL
     nvar = 0
     CALL RGMSG(substr, RGMLI, 'NO VARIABLES LISTED !')
     IF (TRIM(infile) == '') THEN
        CALL RGMSG(substr, RGMLE,'NO INPUT FILE GIVEN !')
     ELSE
        CALL RGMSG(substr, RGMLVM, &
             'SCANNING INPUT FILE '''//TRIM(infile)//''' ...')
        CALL SCAN_NCVAR(kvar, file=TRIM(infile))
        nvar = SIZE(kvar)
        ALLOCATE(ivar(nvar), STAT=status)
        CALL ERRMSG(substr,status,1)
        ALLOCATE(ovar(nvar), STAT=status)
        CALL ERRMSG(substr,status,2)
        ALLOCATE(scl(nvar), STAT=status)
        CALL ERRMSG(substr,status,3)
        ALLOCATE(namatts(nvar), STAT=status)
        CALL ERRMSG(substr,status,4)
        ALLOCATE(RGT(nvar), STAT=status)
        CALL ERRMSG(substr,status,5)
        ALLOCATE(RGTstr(nvar), STAT=status)
        CALL ERRMSG(substr,status,6)
        ALLOCATE(RGrange(nvar,2), STAT=status)
        CALL ERRMSG(substr,status,7)
        !
        ! INIT
        scl(:) = 1.0_dp
        RGrange(:,:) = 0
        RGTrangestr = ''
        !
        DO i=1, nvar
           ivar(i) = TRIM(kvar(i)%name)
           ovar(i) = TRIM(kvar(i)%name)
           ! TRY TO GET RGTstr, RGT from attributes
           RGTstr(i) = '...'
           DO j=1, kvar(i)%natts
              IF (TRIM(kvar(i)%att(j)%name) == 'RG_TYPE') THEN
                 RGTstr(i) = TRIM(string(kvar(i)%att(j)%dat%vc))
              END IF
              IF (TRIM(kvar(1)%att(j)%name) == 'RG_INDEX_RANGE') THEN
                 RGrange(i,1:2) = INT(kvar(1)%att(j)%dat%vi(1:2))
              END IF
           END DO
        END DO
        ! CLEAN UP
        DO i=1, nvar
           CALL INIT_NCVAR(kvar(i))
        END DO
        DEALLOCATE(kvar, STAT=status)
        CALL ERRMSG(substr,status,8)
        NULLIFY(kvar)
     END IF

  ELSE     ! VAR LIST NOT EMPTY

     ! COUNT VARIABLES
     nvar = 1
     len = LEN_TRIM(var)
     DO i=1, len-1        ! last ';' optional !
        IF (var(i:i) == ';') nvar = nvar + 1
     END DO

     ! ALLOCATE SPACE
     ALLOCATE(ivar(nvar), STAT=status)
     CALL ERRMSG(substr,status,9)
     ALLOCATE(ovar(nvar), STAT=status)
     CALL ERRMSG(substr,status,10)
     ALLOCATE(scl(nvar), STAT=status)
     CALL ERRMSG(substr,status,11)
     ALLOCATE(RGT(nvar), STAT=status)
     CALL ERRMSG(substr,status,12)
     ALLOCATE(RGTstr(nvar), STAT=status)
     CALL ERRMSG(substr,status,13)
     ALLOCATE(lvar(nvar), STAT=status)
     CALL ERRMSG(substr,status,14)
     ALLOCATE(namatts(nvar), STAT=status)
     CALL ERRMSG(substr,status,15)
     ALLOCATE(attvar(nvar), STAT=status)
     CALL ERRMSG(substr,status,16)
     ALLOCATE(RGrange(nvar,2), STAT=status)
     CALL ERRMSG(substr,status,17)
     !
     ! INIT
     scl(:)   = 1.0_dp
     RGrange(:,:) = 0
     !
     ! PARSING
     idx1 = 1
     idx2 = INDEX(TRIM(var), ';')
     DO i=1, nvar-1
        idxatt = idx1+INDEX(var(idx1:idx2-1), '|')
        IF (idxatt == idx1) THEN
           idxatt = idx2
           attvar(i) = ' '
        ELSE
           attvar(i) =  var(idxatt:idx2-1)
        END IF
        lvar(i) = var(idx1:idxatt-1)
        idx1 = idx2+1
        idx2 = idx2+INDEX(var(idx1:),';')
     END DO
     IF (idx2 == (idx1-1)) idx2 = len+1    ! last ';' optional !
     idxatt = idx1+INDEX(var(idx1:idx2-1), '|')
     IF (idxatt == idx1) THEN
        idxatt = idx2
        attvar(nvar) = ' '
     ELSE
        attvar(nvar) = var(idxatt:idx2-1)
     END IF
     lvar(nvar)   = var(idx1:idxatt-1)

      ! fill ovar, ivar, scl, RGT, RGTstr
     DO i=1, nvar  ! LOOP OVER ALL VARIABLES
        ! GET indices
        len  = LEN_TRIM(lvar(i))
        idx1 = INDEX(lvar(i),'=')
        IF (idx1 == 0) idx1 = len+1
        idx2 = INDEX(lvar(i),',')
        IF (idx2 == 0) idx2 = len+1
        idx3 = INDEX(lvar(i),':')
        IF (idx3 == 0) idx3 = len+1
        IF (idx1 > idx3) idx1 = len+1 ! filter = in IXF range
        idx4 = INDEX(lvar(i),'|')
        IF (idx4 == 0) idx4 = len+1

        ! CHECK for destination variable name
        idxl = 1
        idxr = MIN(idx1, idx2, idx3, idx4)-1
        ovar(i) = lvar(i)(idxl:idxr)
        ovar(i) = ADJUSTL(TRIM(ovar(i)))

        ! CHECK for source variable name
        idxl = idx1+1
        idxr = MIN(idx2, idx3, idx4)-1
        IF (idxl > len) idxl = 1
        IF (idxr > len) idxr = len
        ivar(i) = lvar(i)(idxl:idxr)
        ivar(i) = ADJUSTL(TRIM(ivar(i)))
        ! CHECK for scaling
        idxl = idx2+1
        idxr = MIN(idx3, idx4) - 1
        IF (idxr < idxl) idxr = idx4-1
        IF (idxl > len) THEN
           scl(i) = 1.0_dp
        ELSE
           READ(lvar(i)(idxl:idxr),*) scl(i)
        END IF
        ! CHECK for RE-GRIDDING TYPE
        idxl = idx3+1
        idxr = idx4-1
        IF (idxr < idxl) idxr = len
        IF (idxl > len) THEN        ! RG_TYPE NOT specified
           RGTstr(i) = '...'
           IF (TRIM(infile) /= '') THEN  ! filename OK
              ALLOCATE(kvar(1), STAT=status)
              CALL ERRMSG(substr,status,16)
              CALL IMPORT_NCVAR(kvar(1), 1, &
                   varname=TRIM(ivar(i)), file=TRIM(infile))
              CALL INIT_NARRAY(kvar(1)%dat)
              DO j=1, kvar(1)%natts
                 IF (TRIM(kvar(1)%att(j)%name) == 'RG_TYPE') THEN
                    RGTstr(i) = TRIM(string(kvar(1)%att(j)%dat%vc))
                 END IF
                 IF (TRIM(kvar(1)%att(j)%name) == 'RG_INDEX_RANGE') THEN
                    RGrange(i,1:2) = INT(kvar(1)%att(j)%dat%vi(1:2))
                 END IF
              END DO
              ! CLEAN UP
              CALL INIT_NCVAR(kvar(1))
              DEALLOCATE(kvar, STAT=status)
              CALL ERRMSG(substr,status,17)
              NULLIFY(kvar)
           ELSE                          ! filename not OK
              CALL RGMSG(substr, RGMLW, &
                   'REGRIDDING TYPE OF VARIABLE '''//TRIM(ivar(i))//'''')
              CALL RGMSG(substr, RGMLWC, &
                   'CANNOT BE DETERMINED FROM FILE, ')
              CALL RGMSG(substr, RGMLWC, &
                   'SINCE FILENAME IS EMPTY/MISSING !')
           ENDIF                         ! filename
        ELSE                      ! RG_TYPE specified
           RGTrangestr = lvar(i)(idxl:idxr)
           len  = LEN_TRIM(RGTrangestr)
           idx1 = INDEX(RGTrangestr,'=')
           IF (idx1 == 0) idx1 = len+1
           RGTstr(i) = RGTrangestr(1:idx1-1)
           RGTstr(i) = ADJUSTL(TRIM(RGTstr(i)))
           ! range will only be set, if IXF regridding requested...
           IF ((RGTstr(i)=='IXF') .OR. (RGTstr(i)=='ixf')) THEN
              idx2 = INDEX(RGTrangestr,':')
              IF (idx2 == 0) idx2 = len+1
              IF (idx1 < len) THEN
                 idxl = idx1+1
                 idxr = idx2-1
                 READ(RGTrangestr(idxl:idxr),*) RGrange(i,1)
              END IF
              IF (idx2 < len) THEN
                 idxl = idx2+1
                 idxr = len
                 READ(RGTrangestr(idxl:idxr),*) RGrange(i,2)
              ELSE
                 RGrange(i,2) = RGrange(i,1)
                 RGrange(i,1) = 1
              END IF
           END IF
           RGTrangestr = ''
        END IF                    ! RG_TYPE (NOT) specified

        ! ANALYSE ATTRIBUTE STRING
        natt = 1
        len = LEN_TRIM(attvar(i))
        IF (len > 0) THEN
           DO ix=1, len-1        ! last ';' optional !
              IF (attvar(i)(ix:ix) == '|') natt = natt + 1
           END DO
           namatts(i)%natts = natt
           ALLOCATE(namatts(i)%atts(natt))

           ALLOCATE(varatt(natt))
           idx1 = 1
           idx2 = INDEX(TRIM(attvar(i)), '|')
           DO ix = 1, natt-1
              varatt(ix) = attvar(i)(idx1:idx2-1)
              idx1 = idx2+1
              idx2 = idx2+INDEX(attvar(i)(idx1:),'|')
           END DO
           IF (idx2 == (idx1-1)) idx2 = len+1    ! last ';' optional !
           varatt(natt) = attvar(i)(idx1:idx2-1)

           DO ix = 1, natt
              str=' '
              len  = LEN_TRIM(varatt(ix))
              idx1 = INDEX(varatt(ix),'=')
              IF (idx1 == 0) CYCLE
              CALL INIT_NCATT(namatts(i)%atts(ix))
              namatts(i)%atts(ix)%name = ADJUSTL(TRIM(varatt(ix)(1:idx1-1)))
              namatts(i)%atts(ix)%ID   = ix
              str = TRIM(ADJUSTL(varatt(ix)(idx1+1:len)))
              namatts(i)%atts(ix)%dat%n = LEN_TRIM(str)
              ALLOCATE(namatts(i)%atts(ix)%dat%vc(namatts(i)%atts(ix)%dat%n))
              CALL ERRMSG(substr,status,18)
              DO jx=1, LEN_TRIM(str)
                 namatts(i)%atts(ix)%dat%vc(jx) = str(jx:jx)
              END DO
              namatts(i)%atts(ix)%dat%type = VTYPE_CHAR
           END DO

           IF (ALLOCATED(varatt)) THEN
              DEALLOCATE(varatt, STAT=status)
              CALL ERRMSG(substr,status,19)
           END IF
        END IF

     END DO  ! LOOP OVER ALL VARIABLES

     ! CLEAN
     IF (ALLOCATED(lvar)) THEN
        DEALLOCATE(lvar, STAT=status)
        CALL ERRMSG(substr,status,20)
     END IF
     IF (ALLOCATED(attvar)) THEN
        DEALLOCATE(attvar, STAT=status)
        CALL ERRMSG(substr,status,21)
     END IF
  END IF   ! VAR LIST EMPTY OR NOT

  DO i=1, nvar
     ! CHECK for RE-GRIDDING TYPE
     SELECT CASE (RGTstr(i))
     CASE ('INT','int')
        RGT(i) = RG_INT
     CASE ('EXT','ext')
        RGT(i) = RG_EXT
     CASE ('IDX','idx')
        RGT(i) = RG_IDX
     CASE ('IXF','ixf')
        RGT(i) = RG_IXF
     CASE DEFAULT
        CALL RGMSG(substr, RGMLW, &
             'UNKNOWN RE-GRIDDING TYPE OF VARAIBLE '''//TRIM(ivar(i))//'''')
        CALL RGMSG(substr, RGMLWC, &
             'USING DEFAULT: INT !')
        RGTstr(i) = 'INT'
        RGT(i) = RG_INT
     END SELECT
  END DO

END SUBROUTINE PARSE_VARSTR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE READ_NAMELIST(iou, status                    &
                         ,gi, gg, outfile               &
                         ,var, i_t, g_t, o_t, lp, lint  &
                         ,nmlstr                        &
                         ,infostr)

  USE messy_main_constants_mem,  ONLY: dp
  USE messy_main_grid_trafo,     ONLY: RGEMPTY, GRIDTRAFOVERS
  USE messy_main_grid_netcdf,    ONLY: GRD_MAXSTRLEN, t_narray, RGMSG &
                                     , RGMLW, RGMLE, RGMLI, RGMLIC  &
                                     , INIT_NARRAY
  USE messy_main_grid,           ONLY: t_geohybgrid, INIT_GEOHYBGRID

  IMPLICIT NONE

  ! I/O
  INTEGER                      , INTENT(IN)  :: iou      ! logical I/O unit
  INTEGER                      , INTENT(OUT) :: status   ! file status
  TYPE (t_geohybgrid)          , INTENT(OUT) :: gi       ! input grid
  TYPE (t_geohybgrid)          , INTENT(OUT) :: gg       ! output grid
  CHARACTER (LEN=GRD_MAXSTRLEN), INTENT(OUT) :: outfile  ! output filename
  CHARACTER (LEN=100*GRD_MAXSTRLEN), INTENT(OUT) :: var  ! variable string
  INTEGER                      , INTENT(OUT) :: i_t(4)   ! input time control
  INTEGER                      , INTENT(OUT) :: g_t(3)   ! grid time control
  INTEGER                      , INTENT(OUT) :: o_t(3)   ! output time control
  LOGICAL                      , INTENT(OUT) :: lp       ! pressure or sigma
  LOGICAL                      , INTENT(OUT) :: lint     ! input time ?
  TYPE (t_narray)              , INTENT(OUT) :: nmlstr   ! namelist as string
  CHARACTER(LEN=STRLEN_ULONG)  , INTENT(OUT) :: infostr  ! info string

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'READ_NAMELIST'
  CHARACTER(GRD_MAXSTRLEN) :: infile            ! input netCDF filename
  CHARACTER(GRD_MAXSTRLEN) :: grdfile           ! grid file (netCDF)
  CHARACTER(GRD_MAXSTRLEN) :: i_latm, i_lonm    ! input latitude/longitude
  CHARACTER(GRD_MAXSTRLEN) :: i_lati, i_loni    ! input latitude/longitude
  CHARACTER(GRD_MAXSTRLEN) :: g_latm, g_lonm    ! output latitude/longitude
  CHARACTER(GRD_MAXSTRLEN) :: g_lati, g_loni    ! output latitude/longitude

  ! rotated coordinates
  CHARACTER(GRD_MAXSTRLEN) :: i_rlatm, i_rlonm    ! input latitude/longitude
  CHARACTER(GRD_MAXSTRLEN) :: i_rlati, i_rloni    ! input latitude/longitude
  CHARACTER(GRD_MAXSTRLEN) :: g_rlatm, g_rlonm    ! output latitude/longitude
  CHARACTER(GRD_MAXSTRLEN) :: g_rlati, g_rloni    ! output latitude/longitude
  ! curvilinear (geographical) coordinates
  CHARACTER(GRD_MAXSTRLEN) :: i_clatm, i_clonm    ! input latitude/longitude
  CHARACTER(GRD_MAXSTRLEN) :: i_clati, i_cloni    ! input latitude/longitude
  CHARACTER(GRD_MAXSTRLEN) :: g_clatm, g_clonm    ! output latitude/longitude
  CHARACTER(GRD_MAXSTRLEN) :: g_clati, g_cloni    ! output latitude/longitude

  ! unstructured (geographical) coordinates
  CHARACTER(GRD_MAXSTRLEN) :: g_ulatm, g_ulonm    ! output latitude/longitude
  CHARACTER(GRD_MAXSTRLEN) :: g_ulati, g_uloni    ! output latitude/longitude
  CHARACTER(GRD_MAXSTRLEN) :: i_hyai, i_hybi      ! input hybrid coeff.
  CHARACTER(GRD_MAXSTRLEN) :: i_hyam, i_hybm      ! input hybrid coeff. (mid)
  CHARACTER(GRD_MAXSTRLEN) :: g_hyai, g_hybi      ! output hybrid coeff.
  CHARACTER(GRD_MAXSTRLEN) :: g_hyam, g_hybm      ! output hybrid coeff. (mid)
  CHARACTER(GRD_MAXSTRLEN) :: i_timei, i_timem    ! input time
  CHARACTER(GRD_MAXSTRLEN) :: g_timei, g_timem    ! output time
  CHARACTER(GRD_MAXSTRLEN) :: i_ps, g_ps          ! input/output surface press.
  CHARACTER(GRD_MAXSTRLEN) :: i_p0, g_p0          ! input/output ref. press.
  CHARACTER(GRD_MAXSTRLEN) :: i_pollon, g_pollon  ! input/output ref. press.
  CHARACTER(GRD_MAXSTRLEN) :: i_pollat, g_pollat  ! input/output ref. press.
  CHARACTER(GRD_MAXSTRLEN) :: i_polgam, g_polgam  ! input/output ref. press.
  LOGICAL                  :: pressure            ! pressure or sigma
  LOGICAL                  :: input_time          ! input time information
  !
  REAL(dp), DIMENSION(2)   :: i_latr, i_lonr
  REAL(dp), DIMENSION(2)   :: i_hyar, i_hybr
  LOGICAL                  :: i_lonc
  LOGICAL                  :: i_rlonc
  LOGICAL                  :: i_clonc
  REAL(dp), DIMENSION(2)   :: g_latr, g_lonr
  REAL(dp), DIMENSION(2)   :: g_hyar, g_hybr
  LOGICAL                  :: g_lonc
  LOGICAL                  :: g_rlonc
  LOGICAL                  :: g_clonc

  ! input/output 3d pressure field: required if grid is defined on height level
  CHARACTER(GRD_MAXSTRLEN) :: i_pressi, g_pressi ! interfaces
  CHARACTER(GRD_MAXSTRLEN) :: i_pressm, g_pressm ! mids
  REAL(dp), DIMENSION(2)   :: i_pressr, g_pressr ! ranges input /output
  CHARACTER(LEN=STRLEN_ULONG) :: info

  ! HELPERS
  CHARACTER(GRD_MAXSTRLEN) :: tstr
  CHARACTER(LEN=8)         :: date
  CHARACTER(LEN=10)        :: time
  CHARACTER(LEN=5)         :: zone

  NAMELIST /regrid/ infile, outfile, grdfile                            &
       , info                                                           &
       , i_latm, i_lati, i_lonm, i_loni                                 &
       , i_hyam, i_hyai, i_hybm, i_hybi                                 &
       , i_timei, i_timem, i_ps, i_p0                                   &
       , g_latm, g_lati, g_lonm, g_loni, g_hyam, g_hyai, g_hybm, g_hybi &
       , g_timei, g_timem, g_ps, g_p0                                   &
       , var, i_t, g_t, o_t, pressure, input_time                       &
       , i_lonc, i_latr, i_lonr, i_hyar, i_hybr                         &
       , g_lonc, g_latr, g_lonr, g_hyar, g_hybr                         &
       , i_rlatm, i_rlati, i_rlonm, i_rloni, g_rlatm                    &
       , g_rlati, g_rlonm, g_rloni                                      &
       , i_clatm, i_clati, i_clonm, i_cloni, g_clatm                    &
       , g_clati, g_clonm, g_cloni                                      &
       , i_rlonc, g_rlonc,i_clonc, g_clonc                              &
       , i_pollon, i_pollat, i_polgam, g_pollon, g_pollat, g_polgam     &
       , g_ulatm, g_ulonm, g_ulati, g_uloni                             &
       , i_pressi, i_pressm, i_pressr, g_pressi, g_pressm,g_pressr

  ! INIT (DEFAULT VALUES)
  infile = ''
  outfile = ''
  grdfile = ''
  info   = ''
  i_latm = ''
  i_lati = ''
  i_lonm = ''
  i_loni = ''

  i_hyam = ''
  i_hyai = ''
  i_hybm = ''
  i_hybi = ''
  i_timei = ''
  i_timem = ''
  i_ps = ''
  i_p0 = ''

  i_rlatm = ''
  i_rlati = ''
  i_rlonm = ''
  i_rloni = ''
  i_clatm = ''
  i_clati = ''
  i_clonm = ''
  i_cloni = ''
  i_pollon = ''
  i_pollat = ''
  i_polgam = ''

  g_latm = ''
  g_lati = ''
  g_lonm = ''
  g_loni = ''
  g_hyam = ''
  g_hyai = ''
  g_hybm = ''
  g_hybi = ''
  g_timei = ''
  g_timem = ''
  g_ps = ''
  g_p0 = ''

  g_rlatm = ''
  g_rlati = ''
  g_rlonm = ''
  g_rloni = ''
  g_clatm = ''
  g_clati = ''
  g_clonm = ''
  g_cloni = ''
  g_pollon = ''
  g_pollat = ''
  g_polgam = ''

  g_ulatm = ''
  g_ulati = ''
  g_ulonm = ''
  g_uloni = ''

  g_pressi = ''
  g_pressm = ''
  i_pressi = ''
  i_pressm = ''

  var = ''
  i_t(:) = (/ 1,1,0,0 /)
  g_t(:) = (/ 1,1,0 /)
  o_t(:) = (/ 1,1,0 /)
  pressure   = .FALSE.  ! DEFAULT IS: SIGMA LEVELS
  input_time = .TRUE.   ! DEFAULT IS: USE INPUT TIME INFO FOR OUTPUT

  i_lonc      = .TRUE.
  i_latr(:)   = RGEMPTY
  i_lonr(:)   = RGEMPTY
  i_hyar(:)   = RGEMPTY
  i_hybr(:)   = RGEMPTY
  i_pressr(:) = RGEMPTY
  g_lonc      = .TRUE.
  g_latr(:)   = RGEMPTY
  g_lonr(:)   = RGEMPTY
  g_hyar(:)   = RGEMPTY
  g_hybr(:)   = RGEMPTY
  g_pressr(:) = RGEMPTY

  i_rlonc    = .TRUE.
  g_rlonc    = .TRUE.
  i_rlonc    = .TRUE.
  g_clonc    = .TRUE.

  CALL INIT_GEOHYBGRID(gi)
  CALL INIT_GEOHYBGRID(gg)

  ! READ NAMELIST
  READ(iou, NML=regrid, IOSTAT=status)
  IF (status /= 0) THEN
     RETURN
  END IF

  ! TIME CONTROL CORRECTIONS
  IF (i_t(2) == 0) THEN
     CALL RGMSG(substr, RGMLW, 'INFILE TIME STEP IS ZERO! RESET TO 1!')
     i_t(2) = 1
  END IF
  IF (g_t(2) == 0) THEN
     CALL RGMSG(substr, RGMLW, 'GRDFILE TIME STEP IS ZERO! RESET TO 1!')
     g_t(2) = 1
  END IF
  IF (o_t(2) == 0) THEN
     CALL RGMSG(substr, RGMLW, 'OUTFILE TIME STEP IS ZERO! RESET TO 1!')
     o_t(2) = 1
  END IF

  ! SAVE NAMELIST IN N-ARRAY
  CALL INIT_NARRAY(nmlstr)
  CALL ADDLINE_CNARRAY(nmlstr, 'IMPORT_GRID VERSION '//TRIM(GRIDTRAFOVERS))
  CALL DATE_AND_TIME(date, time, zone)
  CALL ADDLINE_CNARRAY(nmlstr, 'RG_DATE: '//TRIM(date))
  CALL ADDLINE_CNARRAY(nmlstr, 'RG_TIME: '//TRIM(time)//TRIM(zone))
  CALL ADDLINE_CNARRAY(nmlstr, 'NAMELIST:')
  CALL ADDLINE_CNARRAY(nmlstr, '----------------------------------')
  CALL ADDLINE_CNARRAY(nmlstr, 'infile  = '''//TRIM(infile)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'outfile = '''//TRIM(outfile)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'grdfile = '''//TRIM(grdfile)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_latm  = '''//TRIM(i_latm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_lati  = '''//TRIM(i_lati)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_lonm  = '''//TRIM(i_lonm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_loni  = '''//TRIM(i_loni)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_hyam  = '''//TRIM(i_hyam)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_hyai  = '''//TRIM(i_hyai)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_hybm  = '''//TRIM(i_hybm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_hybi  = '''//TRIM(i_hybi)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_timei = '''//TRIM(i_timei)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_timem = '''//TRIM(i_timem)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_ps    = '''//TRIM(i_ps)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_p0    = '''//TRIM(i_p0)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_pressi = '''//TRIM(i_pressi)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_pressm = '''//TRIM(i_pressm)//''',')

  CALL ADDLINE_CNARRAY(nmlstr, 'i_rlatm  = '''//TRIM(i_rlatm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_rlati  = '''//TRIM(i_rlati)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_rlonm  = '''//TRIM(i_rlonm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_rloni  = '''//TRIM(i_rloni)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_clatm  = '''//TRIM(i_clatm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_clati  = '''//TRIM(i_clati)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_clonm  = '''//TRIM(i_clonm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_cloni  = '''//TRIM(i_cloni)//''',')

  CALL ADDLINE_CNARRAY(nmlstr, 'g_latm  = '''//TRIM(g_latm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_lati  = '''//TRIM(g_lati)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_lonm  = '''//TRIM(g_lonm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_loni  = '''//TRIM(g_loni)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_hyam  = '''//TRIM(g_hyam)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_hyai  = '''//TRIM(g_hyai)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_hybm  = '''//TRIM(g_hybm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_hybi  = '''//TRIM(g_hybi)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_timei = '''//TRIM(g_timei)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_timem = '''//TRIM(g_timem)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_ps    = '''//TRIM(g_ps)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_p0    = '''//TRIM(g_p0)//''',')

  CALL ADDLINE_CNARRAY(nmlstr, 'g_rlatm  = '''//TRIM(g_rlatm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_rlati  = '''//TRIM(g_rlati)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_rlonm  = '''//TRIM(g_rlonm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_rloni  = '''//TRIM(g_rloni)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_clatm  = '''//TRIM(g_clatm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_clati  = '''//TRIM(g_clati)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_clonm  = '''//TRIM(g_clonm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_cloni  = '''//TRIM(g_cloni)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_pressi = '''//TRIM(g_pressi)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_pressm = '''//TRIM(g_pressm)//''',')

  CALL ADDLINE_CNARRAY(nmlstr, 'g_ulatm  = '''//TRIM(g_ulatm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_ulati  = '''//TRIM(g_ulati)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_ulonm  = '''//TRIM(g_ulonm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_uloni  = '''//TRIM(g_uloni)//''',')

  CALL ADDLINE_CNARRAY(nmlstr, 'var     = '''//TRIM(var)//''',')
  !
  WRITE(tstr,'(4(I6,",",1x))') i_t
  CALL ADDLINE_CNARRAY(nmlstr, 'i_t     = '//TRIM(tstr))
  WRITE(tstr,'(3(I6,",",1x))') g_t
  CALL ADDLINE_CNARRAY(nmlstr, 'g_t     = '//TRIM(tstr))
  WRITE(tstr,'(3(I6,",",1x))') o_t
  CALL ADDLINE_CNARRAY(nmlstr, 'o_t     = '//TRIM(tstr))
  !
  IF (pressure) THEN
     CALL ADDLINE_CNARRAY(nmlstr, 'pressure= .TRUE.')
  ELSE
     CALL ADDLINE_CNARRAY(nmlstr, 'pressure= .FALSE.')
  END IF
  IF (input_time) THEN
     CALL ADDLINE_CNARRAY(nmlstr, 'input_time= .TRUE.')
  ELSE
     CALL ADDLINE_CNARRAY(nmlstr, 'input_time= .FALSE.')
  END IF
  !
  IF (i_lonc) THEN
     CALL ADDLINE_CNARRAY(nmlstr, 'i_lonc  = .TRUE.')
  ELSE
     CALL ADDLINE_CNARRAY(nmlstr, 'i_lonc  = .FALSE.')
  ENDIF

  WRITE(tstr,'(e10.4,",",e10.4)') i_lonr(1),i_lonr(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'i_lonr  = '//TRIM(tstr))
  WRITE(tstr,'(e10.4,",",e10.4)') i_latr(1),i_latr(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'i_latr  = '//TRIM(tstr))
  WRITE(tstr,'(e10.4,",",e10.4)') i_hyar(1),i_hyar(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'i_hyar  = '//TRIM(tstr))
  WRITE(tstr,'(e10.4,",",e10.4)') i_hybr(1),i_hybr(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'i_hybr  = '//TRIM(tstr))
  !
  IF (i_rlonc) THEN
     CALL ADDLINE_CNARRAY(nmlstr, 'i_rlonc  = .TRUE.')
  ELSE
     CALL ADDLINE_CNARRAY(nmlstr, 'i_rlonc  = .FALSE.')
  ENDIF

  IF (i_clonc) THEN
     CALL ADDLINE_CNARRAY(nmlstr, 'i_clonc  = .TRUE.')
  ELSE
     CALL ADDLINE_CNARRAY(nmlstr, 'i_clonc  = .FALSE.')
  ENDIF

  CALL ADDLINE_CNARRAY(nmlstr, 'i_pollon = '''//TRIM(i_pollon)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_pollat = '''//TRIM(i_pollat)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_polgam = '''//TRIM(i_polgam)//''',')

  IF (g_lonc) THEN
     CALL ADDLINE_CNARRAY(nmlstr, 'g_lonc  = .TRUE.')
  ELSE
     CALL ADDLINE_CNARRAY(nmlstr, 'g_lonc  = .FALSE.')
  ENDIF
  WRITE(tstr,'(e10.4,",",e10.4)') g_lonr(1),g_lonr(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'g_lonr  = '//TRIM(tstr))
  WRITE(tstr,'(e10.4,",",e10.4)') g_latr(1),g_latr(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'g_latr  = '//TRIM(tstr))
  WRITE(tstr,'(e10.4,",",e10.4)') g_hyar(1),g_hyar(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'g_hyar  = '//TRIM(tstr))
  WRITE(tstr,'(e10.4,",",e10.4)') g_hybr(1),g_hybr(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'g_hybr  = '//TRIM(tstr))
  !
  IF (g_rlonc) THEN
     CALL ADDLINE_CNARRAY(nmlstr, 'g_rlonc  = .TRUE.')
  ELSE
     CALL ADDLINE_CNARRAY(nmlstr, 'g_rlonc  = .FALSE.')
  ENDIF

  IF (g_clonc) THEN
     CALL ADDLINE_CNARRAY(nmlstr, 'g_clonc  = .TRUE.')
  ELSE
     CALL ADDLINE_CNARRAY(nmlstr, 'g_clonc  = .FALSE.')
  ENDIF

  CALL ADDLINE_CNARRAY(nmlstr, 'g_pollon = '''//TRIM(g_pollon)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_pollat = '''//TRIM(g_pollat)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_polgam = '''//TRIM(g_polgam)//''',')

  CALL ADDLINE_CNARRAY(nmlstr, ' info = '''//TRIM(info)//''',')

  CALL ADDLINE_CNARRAY(nmlstr, '----------------------------------')

  ! CHECK FILENAMES
  IF (TRIM(infile) == '') THEN
     CALL RGMSG(substr, RGMLE, 'INFILE REQUIRED IN NAMELIST !')
  END IF
  IF (TRIM(outfile) == '') THEN
     CALL RGMSG(substr, RGMLW, 'OUTFILE MISSING IN NAMELIST ! NO OUTPUT !')
  END IF
  IF (TRIM(grdfile) == '') THEN
     CALL RGMSG(substr, RGMLI, &
          'GRDFILE MISSING IN NAMELIST ! INTERFACE MODUS: ')
     CALL RGMSG(substr, RGMLIC, &
          'g_lati, g_latm, g_loni, g_lonm,')
     CALL RGMSG(substr, RGMLIC, &
          'g_hyai, g_hyam, g_hybi, g_hybm,')
     CALL RGMSG(substr, RGMLIC, &
          'g_timei, g_timem, g_ps, g_p0, g_t')
     CALL RGMSG(substr, RGMLIC, 'IGNORED FROM NAMELIST !')
  END IF

  ! TRANSFER INFORMATIONS
  gi%file       = TRIM(infile)
  gi%latm%name  = TRIM(i_latm)
  gi%lati%name  = TRIM(i_lati)
  gi%lonm%name  = TRIM(i_lonm)
  gi%loni%name  = TRIM(i_loni)
  gi%hyam%name  = TRIM(i_hyam)
  gi%hyai%name  = TRIM(i_hyai)
  gi%hybm%name  = TRIM(i_hybm)
  gi%hybi%name  = TRIM(i_hybi)
  gi%timei%name = TRIM(i_timei)
  gi%timem%name = TRIM(i_timem)
  gi%ps%name    = TRIM(i_ps)
  gi%p0%name    = TRIM(i_p0)
  gi%pressi%name = TRIM(i_pressi)
  gi%pressm%name = TRIM(i_pressm)

  gi%rlatm%name  = TRIM(i_rlatm)
  gi%rlati%name  = TRIM(i_rlati)
  gi%rlonm%name  = TRIM(i_rlonm)
  gi%rloni%name  = TRIM(i_rloni)
  gi%clatm%name  = TRIM(i_clatm)
  gi%clati%name  = TRIM(i_clati)
  gi%clonm%name  = TRIM(i_clonm)
  gi%cloni%name  = TRIM(i_cloni)
  gi%pollon%name = TRIM(i_pollon)
  gi%pollat%name = TRIM(i_pollat)
  gi%polgam%name = TRIM(i_polgam)

  gg%file       = TRIM(grdfile)
  gg%latm%name  = TRIM(g_latm)
  gg%lati%name  = TRIM(g_lati)
  gg%lonm%name  = TRIM(g_lonm)
  gg%loni%name  = TRIM(g_loni)
  gg%hyam%name  = TRIM(g_hyam)
  gg%hyai%name  = TRIM(g_hyai)
  gg%hybm%name  = TRIM(g_hybm)
  gg%hybi%name  = TRIM(g_hybi)
  gg%timei%name = TRIM(g_timei)
  gg%timem%name = TRIM(g_timem)
  gg%ps%name    = TRIM(g_ps)
  gg%p0%name    = TRIM(g_p0)
  gg%pressi%name = TRIM(g_pressi)
  gg%pressm%name = TRIM(g_pressm)

  gg%rlatm%name  = TRIM(g_rlatm)
  gg%rlati%name  = TRIM(g_rlati)
  gg%rlonm%name  = TRIM(g_rlonm)
  gg%rloni%name  = TRIM(g_rloni)
  gg%clatm%name  = TRIM(g_clatm)
  gg%clati%name  = TRIM(g_clati)
  gg%clonm%name  = TRIM(g_clonm)
  gg%cloni%name  = TRIM(g_cloni)
  gg%pollon%name = TRIM(g_pollon)
  gg%pollat%name = TRIM(g_pollat)
  gg%polgam%name = TRIM(g_polgam)
  gg%ulatm%name  = TRIM(g_ulatm)
  gg%ulati%name  = TRIM(g_ulati)
  gg%ulonm%name  = TRIM(g_ulonm)
  gg%uloni%name  = TRIM(g_uloni)

  gi%lonc = i_lonc
  gg%lonc = g_lonc

  gi%rlonc = i_rlonc
  gg%rlonc = g_rlonc
  gi%clonc = i_clonc
  gg%clonc = g_clonc

  lp   = pressure
  lint = input_time

  !  var
  !  i_t
  !  g_t
  !  o_t
  !  outfile

  gi%ranges(1,:) = i_lonr(:)
  gi%ranges(2,:) = i_latr(:)
  gi%ranges(3,:) = i_hyar(:)
  gi%ranges(4,:) = i_hybr(:)
  gg%ranges(1,:) = g_lonr(:)
  gg%ranges(2,:) = g_latr(:)
  gg%ranges(3,:) = g_hyar(:)
  gg%ranges(4,:) = g_hybr(:)

  IF (i_pressr(1) /= RGEMPTY .OR. i_pressr(2) /= RGEMPTY) THEN
     gi%ranges(3,:) = i_pressr(:)
     gi%ranges(4,:) = RGEMPTY
  END IF
  IF (g_pressr(1) /= RGEMPTY .OR. g_pressr(2) /= RGEMPTY) THEN
     gg%ranges(3,:) = g_pressr(:)
     gg%ranges(4,:) = RGEMPTY
  END IF

  infostr = TRIM(info)

END SUBROUTINE READ_NAMELIST
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE ADDLINE_CNARRAY(na, str)

  USE messy_main_grid_netcdf, ONLY: t_narray, INIT_NARRAY, COPY_NARRAY &
                                  , VTYPE_UNDEF, VTYPE_CHAR, RGMSG   &
                                  , RGMLE

  IMPLICIT NONE

  ! I/O
  TYPE (t_narray),  INTENT(INOUT) :: na
  CHARACTER(LEN=*), INTENT(IN)    :: str

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'ADDLINE_CNARRAY'
  INTEGER                     :: n, m ! string length
  INTEGER                     :: i    ! counter
  TYPE (t_narray)             :: nas

  ! LOCAL
  IF (na%type == VTYPE_UNDEF) THEN
     ! NEW
     n = LEN_TRIM(str)
     CALL INIT_NARRAY(na, 1, (/ n+1 /), VTYPE_CHAR)   ! +1 for 'NEWLINE'
     DO i=1, n
        na%vc(i) = str(i:i)
     END DO
     na%vc(n+1) = ACHAR(10)
  ELSE
     IF (na%type == VTYPE_CHAR) THEN
        ! SAVE OLD
        CALL INIT_NARRAY(nas)
        CALL COPY_NARRAY(nas, na)
        n = LEN_TRIM(str)
        m = SIZE(nas%vc)
        CALL INIT_NARRAY(na, 1, (/ m+n+1 /), VTYPE_CHAR ) ! +1 for 'NEWLINE'
        DO i=1, m
           na%vc(i) = nas%vc(i)
        END DO
        DO i=1, n
           na%vc(m+i) = str(i:i)
        END DO
        na%vc(n+m+1) = ACHAR(10)
        CALL INIT_NARRAY(nas)
     ELSE
        CALL RGMSG(substr, RGMLE, 'N-ARRAY MUST BE OF TYPE CHAR !')
     END IF
  END IF

END SUBROUTINE ADDLINE_CNARRAY
! ------------------------------------------------------------------
! ******************************************************************

  SUBROUTINE RGTOOL_READ_NCVAR(status, iou, nmlfile, vname, t, var  &
                              , iipol, ogridid, igridid, SCRIP_ID   &
                              , lrg, lrgx, lrgy, lrgz, lok , oarea  &
                              , convgrid, ldompar, lvarpar          &
                              , infostr, mnc)

    ! PERFORMS ONE RE-GRIDDING (TIME-)STEP (t) FOR ONE VARIABLE
    ! (vname) AND RETURNS THE RESULT OF TYPE NCVAR (var)
    ! AND OPTIONALLY THE GRID (grid)
    ! LOGICAL I/O UNIT  : iou          (for namelist-file)
    ! NAMELIST-FILE     : nmlfile
    ! DIMENSION-SWITCHES: lrgx, lrgy, lrgz
    ! REGRIDDING-SWITCH : lrg
    ! SUCCESS-FLAG      : lok
    !
    ! Author: Patrick Joeckel, MPICH, October 2002
    !         Astrid Kerkweg,  UniMz, 2013 (restrctured + expanded)

    USE messy_main_grid_trafo,       ONLY: GTRF_NONE, GTRF_NRGD, GTRF_SCRP   &
                                         , EXPAND_INPUT_GRID_LON
    USE messy_main_grid_tools,       ONLY: RGTOOL_CONVERT                    &
                                         , RGTOOL_CONVERT_DAT2VAR
    USE messy_main_grid_trafo_nrgd,  ONLY: REGRID_CONTROL
    USE messy_main_grid_trafo_nrgd_base, ONLY: SNREGRID
    USE messy_main_grid_trafo_scrp,  ONLY: t_scrip_data , CALC_SCRIP_WEIGHTS &
                                         , SCRIP_CONTROL, CALC_SCRIPDATA     &
                                         , CONSTRUCT_VERTICAL_AXIS           &
                                         , INTERPOL_GEOHYBGRID_PS
    USE messy_main_grid_netcdf,      ONLY: t_ncvar, COPY_NCVAR, QDEF_NCVAR   &
                                         , INIT_NCVAR, NULL_XTYPE            &
                                         , ERRMSG, RGMLI, RGMSG, t_multinc
    USE messy_main_grid,             ONLY: t_geohybgrid, COPY_GEOHYBGRID     &
                                         , NEW_GEOHYBGRID, LOCATE_GEOHYBGRID &
                                         , INIT_GEOHYBGRID!, PRINT_GEOHYBGRID

    IMPLICIT NONE

    INTRINSIC :: PRESENT, SIZE, TRIM

    ! I/O
    INTEGER, INTENT (OUT)              :: status       ! error status
    INTEGER, INTENT (IN)               :: iou          ! logical I/O unit
    CHARACTER(LEN=*), INTENT(IN)       :: nmlfile      ! namelist-file
    CHARACTER(LEN=*), INTENT(IN)       :: vname        ! variable name
    INTEGER, INTENT (IN)               :: t            ! netCDF time step
    TYPE (t_ncvar), INTENT(OUT)        :: var          ! nc-variable
    INTEGER, INTENT(IN)                :: iipol        ! interpolation method
    INTEGER, INTENT(IN)                :: ogridid      ! output grid id
    INTEGER, INTENT(INOUT)             :: igridid      ! output grid id
    INTEGER, INTENT(INOUT)             :: SCRIP_ID     ! SCRIP DATA ID
    LOGICAL, INTENT(IN),  OPTIONAL     :: lrg          ! regrid really ?
    LOGICAL, INTENT(IN),  OPTIONAL     :: lrgx         ! regrid in x
    LOGICAL, INTENT(IN),  OPTIONAL     :: lrgy         ! regrid in y
    LOGICAL, INTENT(IN),  OPTIONAL     :: lrgz         ! regrid in z
    LOGICAL, INTENT(OUT), OPTIONAL     :: lok          ! OK?
    REAL(dp), DIMENSION(:,:), POINTER, OPTIONAL :: oarea
    TYPE(t_geohybgrid),  OPTIONAL, INTENT(INOUT):: convgrid
    LOGICAL, INTENT(IN), OPTIONAL      :: ldompar ! interpolation runs in domain
                                                  ! parallelisation
    LOGICAL, INTENT(IN), OPTIONAL      :: lvarpar ! parallise over input vars
    CHARACTER(LEN=STRLEN_ULONG), INTENT(OUT), OPTIONAL :: infostr
    TYPE(t_multinc), INTENT(INOUT), OPTIONAL :: mnc ! multi-netcdf descriptor

    ! LOCAL
    ! OUTPUT OF REGRIDDING PROCEDURE
    CHARACTER(LEN=*), PARAMETER           :: substr = 'RGTOOL_READ_NCVAR'
    LOGICAL,  SAVE                        :: lint    ! use input time ?
    INTEGER,        DIMENSION(:), POINTER :: RGT   => NULL() ! regridding type
    TYPE (t_ncvar), DIMENSION(:), POINTER :: rvar  => NULL() ! list of variables
    TYPE (t_ncvar), DIMENSION(:), POINTER :: ovar  => NULL() ! list of variables
    TYPE (t_ncvar), DIMENSION(:), POINTER :: sovar => NULL() ! list of variables
    TYPE (t_geohybgrid)                   :: igrid   ! input grid info
    TYPE (t_geohybgrid)                   :: ogrid   ! output grid info
    TYPE (t_geohybgrid)                   :: intgrid ! intermediate grid info
    INTEGER                               :: tt      ! time step
    LOGICAL                               :: llrg    ! regrid really ?
    INTEGER                               :: pcnt    ! procedure counter
    INTEGER                               :: i !, j, k ! counter
    TYPE(t_scrip_data), POINTER           :: PSD => NULL()
    LOGICAL                               :: llrgx, llrgy, llrgz
    INTEGER                               :: ix
    LOGICAL                               :: lpresax = .FALSE.
    LOGICAL                               :: lwork   = .FALSE.
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: dat  => NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: vdat => NULL()
    INTEGER                               :: nvar, nx, ny, nn, nz
    REAL(dp), DIMENSION(:,:,:), POINTER   :: vax_in  => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER   :: vax_out => NULL()
    INTEGER                               :: xsize, ysize
    LOGICAL                               :: lrot = .FALSE.
    INTEGER                               :: zdim
    REAL(dp), DIMENSION(:), ALLOCATABLE   :: help
    CHARACTER(LEN=STRLEN_ULONG)           :: info_str

    ! INITIALIZE
    igridid = -99
    status  = 3999

    IF (PRESENT(lok)) THEN
       lok = .FALSE.
    END IF
    IF (PRESENT(lrg)) THEN
       llrg = lrg
    ELSE
       llrg = .true.  ! DEFAULT
    END IF

    IF (PRESENT(lrgx)) THEN
       llrgx = lrgx
    ELSE
       llrgx = .true.  ! DEFAULT
    END IF
    IF (PRESENT(lrgy)) THEN
       llrgy = lrgy
    ELSE
       llrgy = .true.  ! DEFAULT
    END IF
    IF (PRESENT(lrgz)) THEN
       llrgz = lrgz
    ELSE
       llrgz = .false.  ! DEFAULT
    END IF
    !
    NULLIFY(rvar)

    IF (PRESENT(convgrid)) CALL INIT_GEOHYBGRID(convgrid)
    !
    ! LOCATE OUTPUT GRID
    CALL INIT_GEOHYBGRID(ogrid)

    IF (ogridid /= -99) THEN
       CALL LOCATE_GEOHYBGRID(status, ogridid, grid=ogrid)
       IF (status /= 0) RETURN
    END IF

    CALL INIT_NCVAR(var)

    !
    ! USE 1st TIME STEP FOR SCANNING, SINCE A NAMELIST FILE
    ! WITH MORE THAN ONE NAMELIST MUST NOT NECESSARILY
    ! CONTAIN TIME AXES OF EQUAL LENGTH ...
    tt = t
    !
    ! COUNT PROCEDURE LEVEL: 0 : START, VARIABLE NOT YET FOUND IN NAMELIST FILE
    !                        1 : VARIABLE FOUND IN NAMELSIT FILE, BUT NOT YET
    !                            PROCESSED
    !                        2 : VARIABLE FOUND IN NAMELIST FILE, PROCESSING
    !                            ALREADY PERFORMED
    pcnt = 0

    ! START REGRIDDING
    ! REGRIDDING CONTROLLED BY INPUT VARIABELS
    status  = 3999
    RG_CTRL = RG_SCAN         ! SCAN-MODE, NO REGRIDDING
    RG_NML  = NML_NEXT        ! READ NEXT NAMELIST FROM FILE
    RG_STATUS = RGSTAT_START  ! STATUS START REGRIDDING

    DO ! endless DO loop (must be terminated with EXIT)

       CALL READ_CONTROL(RG_CTRL, RG_NML, RG_STATUS          &
                           , rvar                            &
                           , igrid                           &
                           , ogrid                           &
                           , RGT                             &
                           , lint                            &
                           , TRIM(nmlfile)                   &
                           , info_str                        &
                           , tc=tt                           &
                           , iounit = iou                    &
                           , ldomainpar=ldompar              &
                           , lvarparallel= lvarpar           &
                           , lpresaxis=lpresax               &
                           , lwork=lwork                     &
                           , mnc=mnc                         &
                           )

       IF (RG_STATUS == RGSTAT_STOP) THEN
          IF (RG_NML==NML_NEXT) THEN  ! still searching for correct namelist
             status = 6000
          ELSE
             status = 0
          END IF
          EXIT
       END IF

       DO i=1, SIZE(rvar)
          IF (TRIM(rvar(i)%name) == TRIM(vname)) THEN
             ! VARIABLE FOUND IN NAMELIST:
             ! SET PROCEDURE LEVEL
             pcnt = pcnt + 1
             EXIT
          END IF
       END DO

       ! VARIABLE NOT YET FOUND IN NAMELIST FILE:
       ! SCAN NEXT NAMELIST IN NAMELIST FILE
       IF (pcnt == 0) THEN
          RG_NML  = NML_NEXT   ! ... TRY NEXT NAMELIST
          RG_CTRL = RG_SCAN    ! ... SCAN FOR VARIABLE

       ELSE
          ! VARIABLE FOUND IN NAMELIST FILE:
          ! SET SWITCHES FOR PROCESSING
          RG_NML  = NML_STAY
          IF (llrg) THEN
             pcnt = pcnt + 1
             RG_CTRL = RG_PROC  ! ... PERFORM REGRIDDING
          ELSE
             RG_CTRL = RG_SCAN  ! ... SCAN AGAIN (= IMPORT RAW DATA)
          END IF

          IF (PRESENT(infostr)) infostr = info_str
       END IF

       ! FILENAME FOUND IN NAMELIST FILE AND PROCESSING FINISHED:
       ! SET SWITCHES TO END LOOP
       IF (RG_CTRL == RG_PROC) THEN
          SELECT CASE(iipol)
          CASE(GTRF_NONE)
             CALL COPY_NCVAR(var, rvar(i))  ! ... RETURN VARIABLE
             IF (PRESENT(convgrid)) CALL COPY_GEOHYBGRID(convgrid, igrid)
          CASE(GTRF_NRGD)
             IF (lwork) THEN
                CALL REGRID_CONTROL(igrid, ogrid, rvar, ovar, RGT, lint   &
                     , lrgx=llrgx, lrgy=llrgy, lrgz=llrgz, lfirsto=.TRUE. &
                     , lpresaxis=lpresax, grid_conv=intgrid)
                CALL COPY_NCVAR(var, ovar(i))  ! ... RETURN VARIABLE
             ENDIF
             IF (PRESENT(convgrid)) THEN
                IF (lwork) THEN
                   CALL COPY_GEOHYBGRID(convgrid, intgrid)
                ELSE
                   CALL COPY_GEOHYBGRID(convgrid, ogrid)
                END IF
             END IF
             CALL INIT_GEOHYBGRID(intgrid)
          CASE(GTRF_SCRP)
             IF (igrid%lonm%xtype /= NULL_XTYPE) THEN
                IF (igrid%lonm%dim(1)%len==1) THEN
                   CALL EXPAND_INPUT_GRID_LON(igrid, rvar)
                END IF
             END IF
             IF (llrgx .OR. llrgy) THEN
                CALL SCRIP_CONTROL(status, SCRIP_ID, igrid, ogrid, RGT, lint &
                     , rvar, sovar, intgrid, llrgz=llrgz, lfirsto=.TRUE.)
                IF (status /= 0) RETURN
             ELSE
                ALLOCATE(sovar(SIZE(rvar)))
                DO ix = 1, SIZE(rvar)
                   CALL COPY_NCVAR(sovar(ix), rvar(ix))
                END DO
                CALL COPY_GEOHYBGRID(intgrid, igrid)
             END IF

             ! TODO vertical interpolation
             ifvert: IF (llrgz) THEN

                ! make grid with vertical input grid, but horizontal
                ! SCRIP-gridded grid

                ! CONVERT VARIABLE into 4D SPACE
                CALL RGMSG(substr, RGMLI, 'CONVERT SOVAR for hori. dims')
                CALL RGTOOL_CONVERT(sovar(1), dat, intgrid,order='xyzn')
                xsize = SIZE(dat,1)
                ysize = SIZE(dat,2)
                DEALLOCATE(dat)
                NULLIFY(dat)

                ! a DEFINITION OF VERTICAL AXIS: IN-GRID
                CALL RGMSG(substr, RGMLI, 'DEFINE VERTICAL InGrid Axis')
                !
                CALL CONSTRUCT_VERTICAL_AXIS(status, xsize, ysize, vax_in &
                     , igrid, lpresax, SCRIP_ID, ogrid, RGT, lint)
                IF (status /= 0) &
                     CALL ERRMSG('CONSTRUCT_VERTICAL_AXIS: ' ,status,1)

                ! a DEFINITION OF VERTICAL AXIS: IN-GRID
                CALL RGMSG(substr, RGMLI, 'DEFINE VERTICAL OutGrid Axis')
                !
                ! CHECK FOR SPECIAL CASE, vertical axis defined with hybrid
                ! coefficients without surface pressure => get surface pressure
                ! from input grid (if available)
                IF (.NOT.(QDEF_NCVAR(ogrid%pressi) &
                     .OR. QDEF_NCVAR(ogrid%pressm))) THEN
                   IF (QDEF_NCVAR(ogrid%hybi)  &
                        .AND.( .NOT. QDEF_NCVAR(ogrid%ps))) THEN
                      IF (QDEF_NCVAR(igrid%ps)) THEN
                         CALL INTERPOL_GEOHYBGRID_PS(status, igrid, ogrid &
                              , SCRIP_ID)
                      ELSE
                         CALL ERRMSG(&
                              'WRONG INPUT FOR VERTICAL AXIS DEFINITION: '&
                              ,42,1)
                      END IF
                   END IF
                END IF

                CALL CONSTRUCT_VERTICAL_AXIS(status, xsize, ysize, vax_out &
                     , ogrid, lpresax)
                IF (status /= 0) &
                     CALL ERRMSG('CONSTRUCT_VERTICAL_AXIS: ' ,status,1)

                ! CHECK if both vertical axis are orientated in the same way
                ! assume axis orientation is equal for all columns, check only
                ! point (1,1)
                IF ( ( (vax_in(1,1,1)-vax_in(1,1,SIZE(vax_in,3))) * &
                     (vax_out(1,1,1)-vax_out(1,1,SIZE(vax_out,3)))) < 0) THEN
                   ! V-AXIS orientation differs
                   ! => AXIS ROTATION in VAX_IN and DAT required
                   lrot = .TRUE.
                ELSE
                   lrot = .FALSE.
                END IF
                numvar: DO nvar =1, SIZE(sovar)
                   ! CONVERT VARIABLE ON INTGRID TO 3D
                   ! NOTE: vertical coordinate is in intgrid 'n'
                   CALL RGTOOL_CONVERT(sovar(nvar), dat, intgrid,order='xyzn')
                   ! ALLOCATE MEMORY FOR VERTICAL REMAPPED FIELD
                   ALLOCATE(vdat(SIZE(dat,1),SIZE(dat,2) &
                        ,SIZE(vax_out,3)-1,SIZE(dat,4)))
                   vdat = 0._dp
                   ! ROTATE VAX-IN and DATA IF REQUIRED
                   IF (lrot) THEN
                      zdim = SIZE(VAX_IN,3)
                      ALLOCATE(help(zdim))
                      DO nx = 1, SIZE(dat,1)
                         DO ny = 1, SIZE(dat,2)
                            help(:) = VAX_IN(nx,ny,:)
                            DO nz = 1, SIZE(help)
                               VAX_IN(nx,ny,nz) = help(zdim-nz+1)
                            END DO
                            DO nn = 1, SIZE(dat,4)
                               help(1:zdim-1) = dat(nx,ny,:,nn)
                               DO nz = 1, zdim-1
                                  dat(nx,ny,nz,nn) = help(zdim-1-nz+1)
                               END DO
                            END DO
                         END DO
                      END DO
                      DEALLOCATE(help)
                   END IF
                   DO nx = 1, SIZE(dat,1)
                      DO ny = 1, SIZE(dat,2)
                         DO nn = 1, SIZE(dat,4)
                            !
                            CALL SNREGRID(status                             &
                                 , vax_in(nx,ny,:), vax_out(nx,ny,:)         &
                                 , dat(nx,ny,:,nn), vdat(nx,ny,:,nn), .FALSE.)
                            IF (status /= 0) &
                                 CALL ERRMSG('SNREGRID ERROR: ' ,status,23)
                         END DO
                      END DO
                   END DO
                   ! convert vdat to var ..
                   CALL RGTOOL_CONVERT_DAT2VAR(var, vdat &
                        , var%name, ogrid, 'xyzn')
                   ! free memory
                   DEALLOCATE(dat)
                   DEALLOCATE(vdat)

                END DO numvar
                DEALLOCATE(vax_in)
                NULLIFY(vax_in)
                DEALLOCATE(vax_out)
                NULLIFY(vax_out)

                IF (PRESENT(convgrid)) CALL COPY_GEOHYBGRID(convgrid, ogrid)

             ELSE
                CALL COPY_NCVAR(var, sovar(1))
                IF (PRESENT(convgrid)) &
                     CALL COPY_GEOHYBGRID(convgrid, intgrid)
             END IF ifvert

             CALL INIT_GEOHYBGRID(intgrid)
          END SELECT
          RG_CTRL = RG_STOP         ! ... TERMINATE ...
                                     ! ... REGRIDDING PROCEDURE CORRECTLY
          IF (PRESENT(lok)) THEN
             lok = .TRUE.                        ! ... REPORT SUCCESS
          END IF
      ELSE IF ( RG_NML==NML_STAY) THEN! IF PROC
          CALL COPY_NCVAR(var, rvar(i))  ! ... RETURN VARIABLE
          ! manipulate input grid mask (important for large data sets)
          !CALL CONSTRUCT_IGRIDMASK(igrid, ogrid)

          CALL NEW_GEOHYBGRID(status, igridid, igrid)
          IF (status /= 0 .AND. status /= 1) RETURN

          IF (iipol == GTRF_SCRP) THEN
             IF (igrid%lonm%xtype /= NULL_XTYPE) THEN
                IF (igrid%lonm%dim(1)%len==1) THEN
                   CALL EXPAND_INPUT_GRID_LON(igrid, rvar, .TRUE.)
                END IF
             END IF
              CALL CALC_SCRIPDATA(status, igrid, ogrid, RGT, SCRIP_ID &
                  , oarea, PSD, label='import')
              IF (status /= 0 ) THEN
                IF (status /= 1 )  RETURN
             ELSE
                ! CALCULATE WEIGHTS
                CALL CALC_SCRIP_WEIGHTS(status, PSD)
                IF (status /=0) RETURN
             END IF
             NULLIFY(PSD)
          ENDIF
          IF (PRESENT(convgrid)) CALL COPY_GEOHYBGRID(convgrid, igrid)
          !
          RG_CTRL = RG_STOP          ! ... TERMINATE ...
                                     ! ... REGRIDDING PROCEDURE CORRECTLY
          IF (PRESENT(lok)) THEN
             lok = .TRUE.                        ! ... REPORT SUCCESS
          END IF
       END IF

       IF (ASSOCIATED(rvar)) THEN
          DO i=1, SIZE(rvar)
             CALL INIT_NCVAR(rvar(i))
          END DO
          DEALLOCATE(rvar)
          NULLIFY(rvar)
       END IF

       IF (ASSOCIATED(ovar)) THEN
          DO i=1, SIZE(ovar)
             CALL INIT_NCVAR(ovar(i))
          END DO
          DEALLOCATE(ovar)
          NULLIFY(ovar)
       ENDIF

      IF (ASSOCIATED(sovar)) THEN
          DO i=1, SIZE(sovar)
             CALL INIT_NCVAR(sovar(i))
          END DO
          DEALLOCATE(sovar)
          NULLIFY(sovar)
       ENDIF

       IF  (ASSOCIATED(RGT)) THEN
          DEALLOCATE(RGT) ; NULLIFY(RGT)
       END IF

    END DO  ! ENDLESS DO-LOOP ...
    ! END REGRIDDING

    ! FREE MEMORY
    CALL INIT_GEOHYBGRID(ogrid)
    CALL INIT_GEOHYBGRID(igrid)

  END SUBROUTINE RGTOOL_READ_NCVAR
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_READ_NCFILE(status, iou, nmlfile, fname, t, var &
                               , iipol, ogridid, igridid, SCRIP_ID  &
                               , lrg, lrgx, lrgy, lrgz, lok, oarea  &
                               , convgrid, ldompar, lvarpar         &
                               , infostr, mnc)

    ! PERFORMS ONE RE-GRIDDING (TIME-)STEP (t) FOR ONE FILE
    ! (fname) AND RETURNS LIST OF TYPE NCVAR (var)
    ! AND OPTIONALLY THE GRID (grid)
    ! LOGICAL I/O UNIT  : iou          (for namelist-file)
    ! NAMELIST-FILE     : nmlfile
    ! DIMENSION-SWITCHES: lrgx, lrgy, lrgz
    ! REGRIDDING-SWITCH : lrg
    ! SUCCESS-FLAG      : lok
    !
    ! IF fname = '' THEN PROCESS 1st NAMELIST ONLY !!!
    !
    ! Author: Patrick Joeckel, MPICH, October 2002
    !         Astrid Kerkweg,  UniMz, 2013 (restrctured + expanded)

    ! NCREGRID INTERFACE
    USE messy_main_grid_trafo,        ONLY: GTRF_NONE, GTRF_NRGD, GTRF_SCRP   &
                                          , EXPAND_INPUT_GRID_LON
    USE messy_main_grid_tools,        ONLY: RGTOOL_CONVERT                    &
                                          , RGTOOL_CONVERT_DAT2VAR
    USE messy_main_grid_trafo_nrgd,      ONLY: REGRID_CONTROL
    USE messy_main_grid_trafo_nrgd_base, ONLY: SNREGRID
    USE messy_main_grid_trafo_scrp,   ONLY: t_scrip_data , CALC_SCRIP_WEIGHTS &
                                          , SCRIP_CONTROL, CALC_SCRIPDATA     &
                                          , CONSTRUCT_VERTICAL_AXIS           &
                                          , INTERPOL_GEOHYBGRID_PS
    USE messy_main_grid_netcdf,       ONLY: t_ncvar, COPY_NCVAR, QDEF_NCVAR   &
                                          , INIT_NCVAR, NULL_XTYPE            &
                                          , ERRMSG, RGMLI, RGMSG, t_multinc
    USE messy_main_grid,              ONLY: t_geohybgrid,   COPY_GEOHYBGRID   &
                                          , NEW_GEOHYBGRID, LOCATE_GEOHYBGRID &
                                          , INIT_GEOHYBGRID

    IMPLICIT NONE

    INTRINSIC :: PRESENT, SIZE, TRIM

    ! I/O
    INTEGER, INTENT (OUT)                 :: status       ! error status
    INTEGER, INTENT (IN)                  :: iou          ! logical I/O unit
    CHARACTER(LEN=*), INTENT(IN)          :: nmlfile      ! namelist-file
    CHARACTER(LEN=*), INTENT(IN)          :: fname        ! filename
    INTEGER, INTENT (IN)                  :: t            ! netCDF time step
    TYPE (t_ncvar), DIMENSION(:), POINTER :: var          ! nc-variables
    INTEGER, INTENT(IN)                   :: iipol        ! interpolation method
    INTEGER, INTENT(INOUT)                :: igridid
    INTEGER, INTENT(IN)                   :: ogridid
    INTEGER, INTENT(INOUT)                :: SCRIP_ID     ! SCRIP DATA ID
    LOGICAL, INTENT(IN),   OPTIONAL       :: lrg          ! regrid really ?
    LOGICAL, INTENT(IN),   OPTIONAL       :: lrgx         ! regrid in x
    LOGICAL, INTENT(IN),   OPTIONAL       :: lrgy         ! regrid in y
    LOGICAL, INTENT(IN),   OPTIONAL       :: lrgz         ! regrid in z
    LOGICAL, INTENT (OUT), OPTIONAL       :: lok          ! OK?
    TYPE(t_geohybgrid), INTENT(INOUT), OPTIONAL :: convgrid
    REAL(dp), DIMENSION(:,:), POINTER, OPTIONAL :: oarea
    LOGICAL, INTENT(IN),   OPTIONAL       :: ldompar ! interpolation runs in
                                                     ! domain parallelisation
    LOGICAL, INTENT(IN),   OPTIONAL       :: lvarpar ! parallise over input vars
    CHARACTER(LEN=STRLEN_ULONG), INTENT(OUT), OPTIONAL :: infostr
    TYPE(t_multinc), INTENT(INOUT), OPTIONAL :: mnc ! multi-netcdf descriptor

    ! LOCAL
    ! OUTPUT OF REGRIDDING PROCEDURE
    CHARACTER(LEN=*), PARAMETER           :: substr = 'RGTOOL_READ_NCFILE'
    LOGICAL,  SAVE                        :: lint    ! use input time ?
    INTEGER, DIMENSION(:),  POINTER       :: RGT => NULL() ! regridding type
    TYPE (t_ncvar), DIMENSION(:), POINTER :: rvar  => NULL() ! list of variables
    TYPE (t_ncvar), DIMENSION(:), POINTER :: ovar  => NULL() ! list of variables
    TYPE (t_ncvar), DIMENSION(:), POINTER :: sovar => NULL() ! list of variables
    TYPE (t_geohybgrid)                   :: igrid   ! output grid info
    TYPE (t_geohybgrid)                   :: ogrid   ! output grid info
    TYPE (t_geohybgrid)                   :: intgrid ! intermediate grid info
    INTEGER                               :: tt      ! time step
    LOGICAL                               :: llrg    ! regrid really ?
    INTEGER                               :: pcnt    ! procedure counter
    INTEGER                               :: i !, j, k ! counter
    TYPE(t_scrip_data), POINTER           :: PSD => NULL()
    LOGICAL                               :: llrgx, llrgy, llrgz
    INTEGER                               :: ix
    LOGICAL                               :: lpresax = .FALSE.
    LOGICAL                               :: lwork   = .FALSE.
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: dat  => NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: vdat => NULL()
    INTEGER                               :: nvar, nx, ny, nn, nz
    REAL(dp), DIMENSION(:,:,:), POINTER   :: vax_in  => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER   :: vax_out => NULL()
    INTEGER                               :: xsize, ysize
    LOGICAL                               :: lrot = .FALSE.
    INTEGER                               :: zdim
    REAL(dp), DIMENSION(:), ALLOCATABLE   :: help
    CHARACTER(LEN=STRLEN_ULONG)           :: info_str

    ! INITIALIZE
    igridid = -99
    status  = 3999

    IF (PRESENT(lok)) THEN
       lok = .FALSE.
    END IF
    IF (PRESENT(lrg)) THEN
       llrg = lrg
    ELSE
       llrg = .true.  ! DEFAULT
    END IF
    IF (PRESENT(lrgx)) THEN
       llrgx = lrgx
    ELSE
       llrgx = .true.  ! DEFAULT
    END IF
    IF (PRESENT(lrgy)) THEN
       llrgy = lrgy
    ELSE
       llrgy = .true.  ! DEFAULT
    END IF
    IF (PRESENT(lrgz)) THEN
       llrgz = lrgz
    ELSE
       llrgz = .false.  ! DEFAULT
    END IF
    !
    NULLIFY(rvar)
    NULLIFY(ovar)
    NULLIFY(sovar)
    !
    IF (PRESENT(convgrid)) CALL INIT_GEOHYBGRID(convgrid)

    ! LOCATE OUTPUT GRID
    CALL INIT_GEOHYBGRID(ogrid)
    IF (ogridid /= -99) THEN
       CALL LOCATE_GEOHYBGRID(status, ogridid, grid=ogrid)
       IF (status /= 0) RETURN
    END IF

    ! USE 1st TIME STEP FOR SCANNING, SINCE A NAMELIST FILE
    ! WITH MORE THAN ONE NAMELIST MUST NOT NECESSARILY
    ! CONTAIN TIME AXES OF EQUAL LENGTH ...
    tt = t
    !
    ! COUNT PROCEDURE LEVEL: 0 : START, FILENAME NOT YET FOUND IN NAMELIST FILE
    !                        1 : FILENAME FOUND IN NAMELSIT FILE, BUT NOT YET
    !                            PROCESSED
    !                        2 : FILENAME FOUND IN NAMELIST FILE, PROCESSING
    !                            ALREADY PERFORMED
    pcnt = 0

    ! START REGRIDDING
    ! REGRIDDING CONTROLLED BY INPUT VARIABELS
    status  = 3999
    RG_CTRL = RG_SCAN         ! SCAN-MODE, NO REGRIDDING
    RG_NML  = NML_NEXT        ! READ NEXT NAMELIST FROM FILE
    RG_STATUS = RGSTAT_START  ! STATUS START REGRIDDING
    !
    DO ! endless DO loop (must be terminated with EXIT)

       CALL READ_CONTROL(RG_CTRL, RG_NML, RG_STATUS        &
                           , rvar                          &
                           , igrid                         &
                           , ogrid                         &
                           , RGT                           &
                           , lint                          &
                           , TRIM(nmlfile)                 &
                           , info_str                      &
                           , tc=tt                         &
                           , iounit = iou                  &
                           , ldomainpar=ldompar            &
                           , lvarparallel= lvarpar         &
                           , lpresaxis=lpresax             &
                           , lwork=lwork                   &
                           , mnc=mnc )

       IF (RG_STATUS == RGSTAT_STOP) THEN
           IF (RG_NML == NML_NEXT) THEN ! still searching for correct namelist
              status = 6000
           ELSE
              status = 0
           END IF
          EXIT
       END IF

       ! FILENAME FOUND IN NAMELIST OR ALREADY PROCESSED
       IF ((TRIM(ogrid%file) == TRIM(fname)).OR.(pcnt == 1) &
            .OR. (TRIM(fname) == '') ) THEN
          ! SET PROCEDURE LEVEL
          pcnt = pcnt + 1
       END IF

       ! FILENAME NOT YET FOUND IN NAMELIST FILE:
       ! SCAN NEXT NAMELIST IN NAMELIST FILE
       IF (pcnt == 0) THEN
          RG_NML  = NML_NEXT   ! ... TRY NEXT NAMELIST
          RG_CTRL = RG_SCAN    ! ... SCAN FOR VARIABLE

       ELSE
          ! FILENAME FOUND IN NAMELIST FILE:
          ! SET SWITCHES FOR PROCESSING
          RG_NML= NML_STAY
          IF (llrg) THEN
             RG_CTRL = RG_PROC  ! ... PERFORM REGRIDDING
          ELSE
             RG_CTRL = RG_SCAN  ! ... COPY RAW DATA (= IMPORT RAW DATA)
          END IF
          IF (PRESENT(infostr)) infostr = info_str
       END IF

       ! FILENAME FOUND IN NAMELIST FILE AND PROCESSING FINISHED:
       ! SET SWITCHES TO END LOOP
       IF (RG_CTRL == RG_PROC) THEN
          SELECT CASE(iipol)
          CASE(GTRF_NONE)
             ALLOCATE(var(SIZE(rvar)), STAT=status)
             CALL ERRMSG(substr,status,1)
             DO i=1, SIZE(rvar)
                CALL COPY_NCVAR(var(i), rvar(i))    ! ... RETURN VARIABLES
             END DO
             IF (PRESENT(convgrid)) CALL COPY_GEOHYBGRID(convgrid, igrid)
          CASE(GTRF_NRGD)
             IF (lwork) THEN
                CALL REGRID_CONTROL(igrid, ogrid, rvar, ovar, RGT, lint  &
                     , lrgx=llrgx, lrgy=llrgy, lrgz=llrgz, lfirsto=.TRUE.&
                     , lpresaxis=lpresax, grid_conv=intgrid)
                ALLOCATE(var(SIZE(ovar)), STAT=status)
                CALL ERRMSG(substr,status,1)
                DO i=1, SIZE(ovar)
                   CALL COPY_NCVAR(var(i), ovar(i))    ! ... RETURN VARIABLES
                END DO
             ENDIF
             IF (PRESENT(convgrid)) THEN
                IF (lwork) THEN
                   CALL COPY_GEOHYBGRID(convgrid, intgrid)
                ELSE
                   CALL COPY_GEOHYBGRID(convgrid, ogrid)
                END IF
             END IF
             CALL INIT_GEOHYBGRID(intgrid)
          CASE(GTRF_SCRP)
             IF (igrid%lonm%xtype /= NULL_XTYPE) THEN
                IF (igrid%lonm%dim(1)%len==1) THEN
                   CALL EXPAND_INPUT_GRID_LON(igrid, rvar)
                END IF
             END IF
             IF (llrgx .OR. llrgy) THEN
                CALL SCRIP_CONTROL(status, SCRIP_ID, igrid, ogrid, RGT, lint &
                     , rvar, sovar,intgrid, llrgz=llrgz, lfirsto=.TRUE.)
                IF (status /= 0) RETURN
             ELSE
                ALLOCATE(sovar(SIZE(rvar)))
                DO ix = 1, SIZE(rvar)
                   CALL COPY_NCVAR(sovar(ix), rvar(ix))
                END DO
                CALL COPY_GEOHYBGRID(intgrid, igrid)
             END IF

             ! vertical interpolation
             IF (llrgz) THEN

                ! make grid with vertical input grid, but horizontal
                ! SCRIP-gridded grid

                ! CONVERT VARIABLE into 4D SPACE
                CALL RGMSG(substr, RGMLI, 'CONVERT SOVAR for hori. dims')
                CALL RGTOOL_CONVERT(sovar(1), dat, intgrid,order='xyzn')
                xsize = SIZE(dat,1)
                ysize = SIZE(dat,2)
                DEALLOCATE(dat)
                NULLIFY(dat)

                ! a DEFINITION OF VERTICAL AXIS: IN-GRID
                CALL RGMSG(substr, RGMLI, 'DEFINE VERTICAL InGrid Axis')
                !
                CALL CONSTRUCT_VERTICAL_AXIS(status, xsize, ysize, vax_in &
                     , igrid, lpresax, SCRIP_ID, ogrid, RGT, lint)
                IF (status /= 0) &
                     CALL ERRMSG('CONSTRUCT_VERTICAL_AXIS: ' ,status,1)

                ! a DEFINITION OF VERTICAL AXIS: IN-GRID
                CALL RGMSG(substr, RGMLI, 'DEFINE VERTICAL OutGrid Axis')
                !
                ! CHECK FOR SPECIAL CASE, vertical axis defined with hybrid
                ! coefficients without surface pressure => get surface pressure
                ! from input grid (if available)
                IF (.NOT.(QDEF_NCVAR(ogrid%pressi) &
                     .OR. QDEF_NCVAR(ogrid%pressm))) THEN
                   IF (QDEF_NCVAR(ogrid%hybi)  &
                        .AND.( .NOT. QDEF_NCVAR(ogrid%ps))) THEN
                      IF (QDEF_NCVAR(igrid%ps)) THEN
                         CALL INTERPOL_GEOHYBGRID_PS(status, igrid, ogrid &
                              , SCRIP_ID)
                      ELSE
                      CALL ERRMSG('WRONG INPUT FOR VERTICAL AXIS DEFINITION: '&
                           ,42,1)
                      END IF
                   END IF
                END IF

                CALL CONSTRUCT_VERTICAL_AXIS(status, xsize, ysize, vax_out &
                     , ogrid, lpresax)
                IF (status /= 0) &
                     CALL ERRMSG('CONSTRUCT_VERTICAL_AXIS: ' ,status,1)
                ! CHECK if both vertical axis are orientated in the same way
                ! assume axis orientation is equal for all columns, check only
                ! point (1,1)
                IF ( ( (vax_in(1,1,1)-vax_in(1,1,SIZE(vax_in,3))) * &
                     (vax_out(1,1,1)-vax_out(1,1,SIZE(vax_out,3)))) < 0) THEN
                   ! V-AXIS orientation differs
                   ! => AXIS ROTATION in VAX_IN and DAT required
                   lrot = .TRUE.
                ELSE
                   lrot = .FALSE.
                END IF
                ALLOCATE(var(SIZE(sovar)))
                numvar: DO nvar =1, SIZE(sovar)
                   ! CONVERT VARIABLE ON INTGRID TO 3D
                   ! NOTE: vertical coordinate is in intgrid 'n'
                   CALL RGTOOL_CONVERT(sovar(nvar), dat, intgrid,order='xyzn')
                   ! ALLOCATE MEMORY FOR VERTICAL REMAPPED FIELD
                   ALLOCATE(vdat(SIZE(dat,1),SIZE(dat,2) &
                        ,SIZE(vax_out,3)-1,SIZE(dat,4)))
                   vdat = 0._dp
                   ! ROTATE VAX-IN and DATA IF REQUIRED
                   IF (lrot) THEN
                      zdim = SIZE(VAX_IN,3)
                      ALLOCATE(help(zdim))
                      DO nx = 1, SIZE(dat,1)
                         DO ny = 1, SIZE(dat,2)
                            help(:) = VAX_IN(nx,ny,:)
                            DO nz = 1, SIZE(help)
                               VAX_IN(nx,ny,nz) = help(zdim-nz+1)
                            END DO
                            DO nn = 1, SIZE(dat,4)
                               help(1:zdim-1) = dat(nx,ny,:,nn)
                               DO nz = 1, zdim-1
                                  dat(nx,ny,nz,nn) = help(zdim-1-nz+1)
                               END DO
                            END DO
                         END DO
                      END DO
                      DEALLOCATE(help)
                   END IF
                   DO nx = 1, SIZE(dat,1)
                      DO ny = 1, SIZE(dat,2)
                         DO nn = 1, SIZE(dat,4)
                            !
                            CALL SNREGRID(status                             &
                                 , vax_in(nx,ny,:), vax_out(nx,ny,:)         &
                                 , dat(nx,ny,:,nn), vdat(nx,ny,:,nn), .FALSE.)
                            IF (status /= 0) &
                                 CALL ERRMSG('SNREGRID ERROR: ' ,status,24)
                         END DO
                      END DO
                   END DO
                   ! convert vdat to var ..
                   CALL RGTOOL_CONVERT_DAT2VAR(var(nvar), vdat &
                        , sovar(nvar)%name, ogrid, 'xyzn')
                   ! free memory
                   DEALLOCATE(dat)
                   DEALLOCATE(vdat)

                END DO numvar
                DEALLOCATE(vax_in)
                NULLIFY(vax_in)
                DEALLOCATE(vax_out)
                NULLIFY(vax_out)

                IF (PRESENT(convgrid)) CALL COPY_GEOHYBGRID(convgrid, ogrid)

             ELSE

                ALLOCATE(var(SIZE(sovar)), STAT=status)
                CALL ERRMSG(substr,status,1)
                DO i=1, SIZE(sovar)
                   CALL COPY_NCVAR(var(i), sovar(i))   ! ... RETURN VARIABLES
                END DO
                IF (PRESENT(convgrid)) &
                     CALL COPY_GEOHYBGRID(convgrid, intgrid)

             ENDIF
             CALL INIT_GEOHYBGRID(intgrid)
          END SELECT

          RG_CTRL = RG_STOP              ! ... TERMINATE ...
                                         ! ... REGRIDDING PROCEDURE CORRECTLY
          IF (PRESENT(lok)) THEN
             lok = .TRUE.                        ! ... REPORT SUCCESS
          END IF
       ELSE IF (RG_NML == NML_STAY) THEN
          ALLOCATE(var(SIZE(rvar)), STAT=status)
          CALL ERRMSG(substr,status,1)
          DO i=1, SIZE(rvar)
             CALL COPY_NCVAR(var(i), rvar(i))    ! ... RETURN VARIABLES
          END DO

          ! manipulate input grid mask (important for large data sets)
          !CALL CONSTRUCT_IGRIDMASK(igrid, ogrid)

          CALL NEW_GEOHYBGRID(status, igridid, igrid)
          IF (status /= 0 .AND. status /= 1) RETURN
          IF (iipol == GTRF_SCRP) THEN
             IF (igrid%lonm%xtype /= NULL_XTYPE) THEN
                IF (igrid%lonm%dim(1)%len==1) THEN
                   CALL EXPAND_INPUT_GRID_LON(igrid, rvar, .TRUE.)
                END IF
             END IF
             CALL CALC_SCRIPDATA(status, igrid, ogrid, RGT, SCRIP_ID &
                  , oarea, PSD, label='import')
             IF (status /= 0 ) THEN
                ! CALCULATE WEIGHTS
                IF (status /= 01 )  RETURN
             ELSE
                CALL CALC_SCRIP_WEIGHTS(status, PSD)
                IF (status /= 0) RETURN
             END IF
             NULLIFY(PSD)
          ENDIF
          IF (PRESENT(convgrid)) CALL COPY_GEOHYBGRID(convgrid, igrid)

          RG_CTRL = RG_STOP              ! ... TERMINATE ...
                                         ! ... REGRIDDING PROCEDURE CORRECTLY
          IF (PRESENT(lok)) THEN
             lok = .TRUE.                        ! ... REPORT SUCCESS
          END IF
       END IF

       IF (ASSOCIATED(rvar)) THEN
          DO i=1, SIZE(rvar)
             CALL INIT_NCVAR(rvar(i))
          END DO
          DEALLOCATE(rvar)
          NULLIFY(rvar)
       END IF

       IF (ASSOCIATED(ovar)) THEN
          DO i=1, SIZE(ovar)
             CALL INIT_NCVAR(ovar(i))
          END DO
          DEALLOCATE(ovar)
          NULLIFY(ovar)
       ENDIF

      IF (ASSOCIATED(sovar)) THEN
          DO i=1, SIZE(sovar)
             CALL INIT_NCVAR(sovar(i))
          END DO
          DEALLOCATE(sovar)
          NULLIFY(sovar)
       ENDIF

       IF (ASSOCIATED(RGT)) THEN
          DEALLOCATE(RGT)
          NULLIFY(RGT)
       END IF

    END DO  ! ENDLESS DO-LOOP ...
    ! END REGRIDDING

    ! FREE MEMORY
    CALL INIT_GEOHYBGRID(ogrid)
    CALL INIT_GEOHYBGRID(igrid)

  END SUBROUTINE RGTOOL_READ_NCFILE

!!$  SUBROUTINE CONSTRUCT_IGRIDMASK(xgi, xgo)
!!$
!!$    USE messy_main_grid_netcdf,  ONLY: t_narray, QDEF_NCVAR, INIT_NARRAY     &
!!$                                     , DOUBLE_NARRAY, COPY_NARRAY, ADD_NCATT &
!!$                                     , NULL_DIMID, NULL_VARID, NULL_XTYPE    &
!!$                                     , VTYPE_INT, errmsg, NF90_INT
!!$    USE messy_main_grid,         ONLY: COPY_GEOHYBGRID, INIT_GEOHYBGRID  &
!!$                                     , t_geohybgrid
!!$
!!$    IMPLICIT NONE
!!$
!!$    TYPE(t_geohybgrid), INTENT(INOUT)  :: xgi
!!$    TYPE(t_geohybgrid), INTENT(IN)     :: xgo
!!$
!!$    ! LOCAL
!!$    TYPE(t_geohybgrid)       :: xgg
!!$    TYPE(t_narray)           :: llonm
!!$    TYPE(t_narray)           :: llatm
!!$    TYPE(t_narray)           :: lloni
!!$    TYPE(t_narray)           :: llati
!!$    REAL(dp)                 :: lonmin, lonmax, latmin, latmax
!!$    INTEGER                  :: dimlenlon, dimlenlat
!!$    INTEGER                  :: ix, jx, dim1, iy
!!$    INTEGER                  :: status
!!$    LOGICAL                  :: l_mids        = .FALSE.
!!$    LOGICAL                  :: l_interfaces  = .FALSE.
!!$    LOGICAL                  :: l_curvilinear = .FALSE.
!!$    REAL(dp)                 :: lon_ranges(2) = (/0._dp,360._dp/)
!!$    TYPE(t_narray)           :: lloni_out
!!$    TYPE(t_narray)           :: llati_out
!!$    LOGICAL                  :: loutirr = .FALSE.
!!$    LOGICAL                  :: lcalc   = .FALSE.
!!$    REAL(dp), PARAMETER      :: dx1     = 2._dp
!!$    REAL(dp), PARAMETER      :: dx2     = 4.0_dp
!!$    REAL(dp)                 :: dx
!!$    REAL(dp), PARAMETER      :: dgap    = 15._dp
!!$    REAL(dp), PARAMETER      :: dlonlat = 4._dp
!!$    REAL(dp)                 :: fac     = 1._dp
!!$
!!$    IF (QDEF_NCVAR(xgi%imask)) THEN
!!$       write (*,*) 'CONSTRUCT_IMASK: imask exists already'
!!$       RETURN
!!$    END IF
!!$
!!$    CALL COPY_GEOHYBGRID(xgg, xgo)
!!$
!!$    dx = dlonlat
!!$
!!$    ! only for unstructured out grids:
!!$    lon_ranges = (/0._dp,360._dp/)
!!$    IF (QDEF_NCVAR(xgg%uloni)) THEN
!!$       CALL DOUBLE_NARRAY(xgg%uloni%dat)
!!$       CALL DOUBLE_NARRAY(xgg%ulati%dat)
!!$       lonmin = MINVAL(xgg%uloni%dat%vd) - dx
!!$       lonmax = MAXVAL(xgg%uloni%dat%vd) + dx
!!$       latmin = MINVAL(xgg%ulati%dat%vd) - dx
!!$       latmax = MAXVAL(xgg%ulati%dat%vd) + dx
!!$       CALL COPY_NARRAY(lloni_out, xgg%uloni%dat)
!!$       CALL DOUBLE_NARRAY(lloni_out)
!!$       CALL COPY_NARRAY(llati_out, xgg%ulati%dat)
!!$       CALL DOUBLE_NARRAY(llati_out)
!!$       loutirr = .TRUE.
!!$       IF (MINVAL(xgg%uloni%dat%vd) < 0) lon_ranges = (/-180._dp , 180._dp/)
!!$    ELSE IF (QDEF_NCVAR(xgg%loni)) THEN
!!$       CALL DOUBLE_NARRAY(xgg%loni%dat)
!!$       CALL DOUBLE_NARRAY(xgg%lati%dat)
!!$       lonmin = MINVAL(xgg%loni%dat%vd) - dx
!!$       lonmax = MAXVAL(xgg%loni%dat%vd) + dx
!!$       latmin = MINVAL(xgg%lati%dat%vd) - dx
!!$       latmax = MAXVAL(xgg%lati%dat%vd) + dx
!!$       loutirr = .FALSE.
!!$       IF (MINVAL(xgg%loni%dat%vd) < 0) lon_ranges = (/-180._dp , 180._dp/)
!!$    ELSE IF (QDEF_NCVAR(xgg%cloni)) THEN
!!$       CALL DOUBLE_NARRAY(xgg%cloni%dat)
!!$       CALL DOUBLE_NARRAY(xgg%clati%dat)
!!$       lonmin = MINVAL(xgg%cloni%dat%vd) - dx
!!$       lonmax = MAXVAL(xgg%cloni%dat%vd) + dx
!!$       latmin = MINVAL(xgg%clati%dat%vd) - dx
!!$       latmax = MAXVAL(xgg%clati%dat%vd) + dx
!!$       CALL COPY_NARRAY(lloni_out, xgg%cloni%dat)
!!$       CALL DOUBLE_NARRAY(lloni_out)
!!$       CALL COPY_NARRAY(llati_out, xgg%clati%dat)
!!$       CALL DOUBLE_NARRAY(llati_out)
!!$       loutirr = .TRUE.
!!$       IF (MINVAL(xgg%cloni%dat%vd) < 0) lon_ranges = (/-180._dp , 180._dp/)
!!$    ELSE
!!$       write (*,*) 'CONSTRUCT_IMASK outgrid not defined sufficiently'
!!$       CALL clean_locals
!!$       RETURN
!!$    END IF
!!$
!!$    CALL INIT_NARRAY(llonm)
!!$    CALL INIT_NARRAY(lloni)
!!$    CALL INIT_NARRAY(llatm)
!!$    CALL INIT_NARRAY(llati)
!!$
!!$    l_mids        = .FALSE.
!!$    l_interfaces  = .FALSE.
!!$    l_curvilinear = .FALSE.
!!$    IF (xgi%lonm%xtype /= NULL_XTYPE) THEN
!!$       ! FORCE xgi%lonm%dat%vd to be associated from now on
!!$       l_mids = .TRUE.
!!$       CALL COPY_NARRAY(llonm, xgi%lonm%dat)
!!$       CALL COPY_NARRAY(llatm, xgi%latm%dat)
!!$    ENDIF
!!$    IF (xgi%loni%xtype /= NULL_XTYPE) THEN
!!$       l_interfaces  = .TRUE.
!!$       CALL COPY_NARRAY(lloni, xgi%loni%dat)
!!$       CALL COPY_NARRAY(llati, xgi%lati%dat)
!!$    ENDIF
!!$
!!$    ! just search curvilinearfields if local llon and llat not yet defined
!!$    IF (.NOT. l_mids .AND. .NOT. l_interfaces) THEN
!!$       IF  (xgi%clonm%xtype /= NULL_XTYPE) THEN
!!$          l_mids        = .TRUE.
!!$          l_curvilinear = .TRUE.
!!$          CALL COPY_NARRAY(llonm, xgi%clonm%dat)
!!$          CALL COPY_NARRAY(llatm, xgi%clatm%dat)
!!$       END IF
!!$       IF (xgi%cloni%xtype /= NULL_XTYPE) THEN
!!$          l_interfaces  = .TRUE.
!!$          l_curvilinear = .TRUE.
!!$          CALL COPY_NARRAY(lloni, xgi%cloni%dat)
!!$          CALL COPY_NARRAY(llati, xgi%clati%dat)
!!$       ENDIF
!!$    END IF
!!$    IF  (.NOT. l_mids .AND. .NOT. l_interfaces) THEN
!!$       CALL clean_locals
!!$       RETURN
!!$    END IF
!!$
!!$    IF (l_mids) THEN
!!$       CALL DOUBLE_NARRAY(llonm)
!!$       CALL DOUBLE_NARRAY(llatm)
!!$       ! SMALL INPUT FIELD DO NOT NEED A MASK
!!$       IF (SIZE(llonm%vd) < 150 ) THEN
!!$          CALL clean_locals
!!$          RETURN
!!$       END IF
!!$    END IF
!!$    IF (l_interfaces) THEN
!!$       CALL DOUBLE_NARRAY(lloni)
!!$       CALL DOUBLE_NARRAY(llati)
!!$       ! SMALL INPUT FIELD DO NOT NEED A MASK
!!$       IF (SIZE(lloni%vd) < 150 ) THEN
!!$          CALL clean_locals
!!$          RETURN
!!$       END IF
!!$    END IF
!!$
!!$    IF (l_interfaces) THEN
!!$       IF (.NOT. l_curvilinear) THEN
!!$          dimlenlon = lloni%dim(1) - 1
!!$          dimlenlat = llati%dim(1) - 1
!!$       ELSE
!!$          dimlenlon = lloni%dim(1) - 1
!!$          dimlenlat = llati%dim(2) - 1
!!$       END IF
!!$       IF (MINVAL(lloni%vd)<= lon_ranges(1)) THEN
!!$          DO ix = 1,SIZE(lloni%vd)
!!$             IF (lloni%vd(ix) < 0._dp) THEN
!!$                lloni%vd(ix) = lloni%vd(ix) + 360._dp
!!$             END IF
!!$          END DO
!!$       ELSE IF (MAXVAL(lloni%vd) > lon_ranges(2)) THEN
!!$          DO ix = 1,SIZE(llonm%vd)
!!$             IF (lloni%vd(ix) > 180._dp) THEN
!!$                lloni%vd(ix) = lloni%vd(ix) - 360._dp
!!$             END IF
!!$          END DO
!!$       END IF
!!$    ELSE ! lmids
!!$       IF (.NOT. l_curvilinear) THEN
!!$          dimlenlon = llonm%dim(1)
!!$          dimlenlat = llatm%dim(1)
!!$       ELSE
!!$          dimlenlon = llonm%dim(1)
!!$          dimlenlat = llatm%dim(2)
!!$       END IF
!!$       IF (MINVAL(llonm%vd)<= lon_ranges(1)) THEN
!!$          DO ix = 1,SIZE(llonm%vd)
!!$             IF (llonm%vd(ix) < 0._dp) THEN
!!$                llonm%vd(ix) = llonm%vd(ix) + 360._dp
!!$             END IF
!!$          END DO
!!$       ELSE IF (MAXVAL(llonm%vd) > lon_ranges(2)) THEN
!!$          DO ix = 1,SIZE(llonm%vd)
!!$             IF (llonm%vd(ix) >  180._dp) THEN
!!$                llonm%vd(ix) = llonm%vd(ix) - 360._dp
!!$             END IF
!!$          END DO
!!$       END IF
!!$    END IF
!!$
!!$    ! DEFINE MASK of defined points
!!$    xgi%imask%name = 'mask'
!!$    xgi%imask%id    = NULL_VARID
!!$    xgi%imask%xtype = NF90_INT
!!$    ! ... dimensions
!!$    xgi%imask%ndims = 2
!!$    ALLOCATE(xgi%imask%dim(xgi%imask%ndims), STAT=status)
!!$    CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
!!$    xgi%imask%dim(1)%name = 'lon'
!!$    xgi%imask%dim(1)%id    = NULL_DIMID
!!$    xgi%imask%dim(1)%len   =  dimlenlon
!!$    xgi%imask%dim(1)%fuid  = .false.
!!$    xgi%imask%dim(1)%varid = NULL_VARID
!!$    xgi%imask%dim(2)%name  = 'lat'
!!$    xgi%imask%dim(2)%id    = NULL_DIMID
!!$    xgi%imask%dim(2)%len   = dimlenlat
!!$    xgi%imask%dim(2)%fuid  = .false.
!!$    xgi%imask%dim(2)%varid = NULL_VARID
!!$
!!$    ! .. ATTRIBUTES
!!$    CALL ADD_NCATT(xgi%imask, 'long_name', vs='mask of defined points')
!!$
!!$    ! ... data
!!$    CALL INIT_NARRAY(xgi%imask%dat, xgi%imask%ndims &
!!$         , (/xgi%imask%dim(1)%len,xgi%imask%dim(2)%len/) &
!!$         ,VTYPE_INT)
!!$
!!$    xgi%imask%dat%vi = 1
!!$
!!$    ! the large dx are not required for small grids
!!$    ! therefore the gap and the distance are reduced by a factor
!!$    fac = 1._dp
!!$    IF (dimlenlon * dimlenlat > 300000._dp) fac = 0.5_dp
!!$
!!$    IF (.NOT. l_curvilinear) THEN
!!$       IF (l_interfaces) THEN
!!$          dim1 = SIZE(lloni%vd) - 1
!!$          DO ix = 1, SIZE(lloni%vd) - 1
!!$             DO jx = 1, SIZE(llati%vd) - 1
!!$                IF (llati%vd(jx) < -90._dp + dgap * fac) CYCLE
!!$                IF (llati%vd(jx) >  90._dp - dgap * fac) CYCLE
!!$                IF (lloni%vd(ix)    < lonmin  .OR. &
!!$                     llati%vd(jx)   < latmin  .OR. &
!!$                     lloni%vd(ix+1) > lonmax  .OR. &
!!$                     llati%vd(jx+1) > latmax)          THEN
!!$                   xgi%imask%dat%vi((jx-1)*dim1 + ix ) = 0
!!$                   CYCLE
!!$                END IF
!!$                IF (loutirr) THEN
!!$                   IF (ABS(llati%vd(jx))  > 60._dp) THEN
!!$                      dx = dx2 * fac
!!$                   ELSE
!!$                      dx = dx1 * fac
!!$                   END IF
!!$                   lcalc = .FALSE.
!!$                   DO iy = 1, SIZE(lloni_out%vd)
!!$                      IF ( ABS(lloni%vd(ix) - lloni_out%vd(iy)) < dx .AND. &
!!$                           ABS(llatm%vd(jx) - llati_out%vd(iy)) < dx) THEN
!!$                         lcalc = .TRUE.
!!$                         EXIT
!!$                      END IF
!!$                   END DO
!!$                   IF (.NOT. lcalc) xgi%imask%dat%vi((jx-1)*dim1 + ix ) = 0
!!$                END IF
!!$             END DO
!!$          END DO
!!$
!!$       ELSE  ! lmids
!!$          dim1 = SIZE(llonm%vd)
!!$          DO ix = 1, SIZE(llonm%vd)
!!$             DO jx = 1, SIZE(llatm%vd)
!!$                IF (llatm%vd(jx) < -90._dp + dgap * fac) CYCLE
!!$                IF (llatm%vd(jx) >  90._dp - dgap * fac) CYCLE
!!$                IF (llonm%vd(ix)  < lonmin  .OR. &
!!$                     llatm%vd(jx) < latmin  .OR. &
!!$                     llonm%vd(ix) > lonmax  .OR. &
!!$                     llatm%vd(jx) > latmax)  THEN
!!$                    xgi%imask%dat%vi((jx-1)*dim1 + ix ) = 0
!!$                   CYCLE
!!$                END IF
!!$                IF (loutirr) THEN
!!$                   lcalc = .FALSE.
!!$                   IF (ABS(llatm%vd(jx))  > 60._dp) THEN
!!$                      dx = dx2 * fac
!!$                   ELSE
!!$                      dx = dx1 * fac
!!$                   END IF
!!$                    DO iy = 1, SIZE(lloni_out%vd)
!!$                      IF ( ABS(llonm%vd(ix) - lloni_out%vd(iy)) < dx .AND. &
!!$                           ABS(llatm%vd(jx) - llati_out%vd(iy))  < dx) THEN
!!$                         lcalc = .TRUE.
!!$                         EXIT
!!$                      END IF
!!$                   END DO
!!$                   IF (.NOT. lcalc) THEN
!!$                      xgi%imask%dat%vi((jx-1)*dim1 + ix ) = 0
!!$                   END IF
!!$                END IF
!!$             END DO
!!$          END DO
!!$       ENDIF
!!$    ELSE ! curvilinear
!!$       IF (l_mids) THEN
!!$          DO ix = 1, SIZE(llonm%vd)
!!$             IF (llatm%vd(ix) < -90._dp + dgap * fac) CYCLE
!!$             IF (llatm%vd(ix) >  90._dp - dgap * fac) CYCLE
!!$             IF (llonm%vd(ix)  < lonmin  .OR. &
!!$                  llatm%vd(ix) < latmin  .OR. &
!!$                  llonm%vd(ix) > lonmax  .OR. &
!!$                  llatm%vd(ix) > latmax)  THEN
!!$                xgi%imask%dat%vi(ix) = 0
!!$                CYCLE
!!$             END IF
!!$
!!$             IF (loutirr) THEN
!!$                IF (ABS(llatm%vd(ix))  > 60._dp) THEN
!!$                   dx = dx2 * fac
!!$                ELSE
!!$                   dx = dx1 * fac
!!$                END IF
!!$                lcalc = .FALSE.
!!$                DO iy = 1, SIZE(lloni_out%vd)
!!$                   IF ( ABS(llonm%vd(ix) - lloni_out%vd(iy)) < dx .AND. &
!!$                        ABS(llatm%vd(ix) - llati_out%vd(iy))  < dx) THEN
!!$                      lcalc = .TRUE.
!!$                      EXIT
!!$                   END IF
!!$                END DO
!!$                IF (.NOT. lcalc) xgi%imask%dat%vi(ix ) = 0
!!$             END IF
!!$          END DO
!!$       ELSE ! interfaces
!!$          DO ix = 1, SIZE(lloni%vd)
!!$             IF (llati%vd(ix) < -90._dp + dgap * fac) CYCLE
!!$             IF (llati%vd(ix) >  90._dp - dgap * fac) CYCLE
!!$             IF (lloni%vd(ix) < lonmin  .OR. &
!!$                  llati%vd(ix) < latmin  .OR. &
!!$                  lloni%vd(ix) > lonmax  .OR. &
!!$                  llati%vd(ix) > latmax)   THEN
!!$                xgi%imask%dat%vi(ix) = 0
!!$                CYCLE
!!$             END IF
!!$
!!$             IF (loutirr) THEN
!!$                lcalc = .FALSE.
!!$                IF (ABS(llati%vd(ix))  > 60._dp) THEN
!!$                   dx = dx2 * fac
!!$                ELSE
!!$                   dx = dx1 * fac
!!$                END IF
!!$                DO iy = 1, SIZE(lloni_out%vd)
!!$                   IF ( ABS(lloni%vd(ix) - lloni_out%vd(iy)) < dx .AND. &
!!$                        ABS(llati%vd(ix) - llati_out%vd(iy))  < dx) THEN
!!$                      lcalc = .TRUE.
!!$                      EXIT
!!$                   END IF
!!$                END DO
!!$                IF (.NOT. lcalc) xgi%imask%dat%vi(ix ) = 0
!!$             END IF
!!$          END DO
!!$       END IF
!!$    END IF
!!$
!!$    CALL clean_locals
!!$
!!$    CONTAINS
!!$
!!$      SUBROUTINE clean_locals
!!$
!!$        CALL INIT_GEOHYBGRID(xgg)
!!$        CALL INIT_NARRAY(llonm)
!!$        CALL INIT_NARRAY(lloni)
!!$        CALL INIT_NARRAY(llatm)
!!$        CALL INIT_NARRAY(llati)
!!$        CALL INIT_NARRAY(lloni_out)
!!$        CALL INIT_NARRAY(llati_out)
!!$
!!$      END SUBROUTINE clean_locals
!!$
!!$  END SUBROUTINE construct_igridmask

! ******************************************************************
! ------------------------------------------------------------------
END MODULE MESSY_MAIN_IMPORT_GRID
! ------------------------------------------------------------------
! ******************************************************************
