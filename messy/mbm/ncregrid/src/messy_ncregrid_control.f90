! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_NCREGRID_CONTROL
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, June 2002
! ******************************************************************

  USE MESSY_NCREGRID_BASE
  USE MESSY_NCREGRID_NETCDF
  USE MESSY_NCREGRID_GEOHYB
!  USE MESSY_NCREGRID_DIAG

  IMPLICIT NONE

  INTRINSIC :: DATE_AND_TIME, TRIM, ACHAR, LEN_TRIM, SIZE      &
             , ALLOCATED, INDEX, MIN, PRESENT, ABS, ASSOCIATED &
             , TINY

  PRIVATE   :: DATE_AND_TIME, TRIM, ACHAR, LEN_TRIM, SIZE      &
             , ALLOCATED, INDEX, MIN, PRESENT, ABS, ASSOCIATED &
             , TINY

#ifndef PNCREGRID
  CHARACTER(*), PARAMETER :: NCREGRIDVERS = '1.5b'
#else
  CHARACTER(*), PARAMETER :: NCREGRIDVERS = '1.5b-mpi'
#endif
  INTEGER, SAVE      :: RG_CTRL, RG_NML, RG_STATUS
  INTEGER, PARAMETER :: RG_SCAN = 1
  INTEGER, PARAMETER :: RG_PROC = 2
  INTEGER, PARAMETER :: RG_STOP = 3
  INTEGER, PARAMETER :: NML_NEXT = 4
  INTEGER, PARAMETER :: NML_STAY = 5
  INTEGER, PARAMETER :: RGSTAT_START = 6
  INTEGER, PARAMETER :: RGSTAT_CONT = 7
  INTEGER, PARAMETER :: RGSTAT_STOP = 8
#ifdef PNCREGRID
  INTEGER, PARAMETER, PUBLIC :: RGSTAT_NULL = -1
  INTEGER, PARAMETER, PUBLIC :: RGSTAT_CYCLE = 9

  TYPE t_mpi_def
     INTEGER :: rank
     INTEGER :: nproc
     INTEGER :: comm
  END TYPE t_mpi_def
  PUBLIC :: t_mpi_def

  INTEGER, PRIVATE, SAVE      :: me, nproc, comm
  INTEGER, PRIVATE, PARAMETER :: MAX_working_PEs=16
  INTEGER, PRIVATE, SAVE      :: nproc_work
  LOGICAL, PRIVATE, SAVE      :: i_am_worker
  INTEGER, PRIVATE, DIMENSION(:), ALLOCATABLE, SAVE :: pe_list
  INTEGER, PRIVATE, SAVE      :: work_loop = 0

  INTERFACE REGRID_CONTROL
     MODULE PROCEDURE REGRID_CONTROL
  END INTERFACE
#endif

CONTAINS

! ------------------------------------------------------------------
SUBROUTINE REGRID_CONTROL(GCTRL, GNML, GSTAT              &
                          ,nmlfile                        &
                          ,iounit                         &
                          ,lrgx, lrgy, lrgz               &
                          ,tmin, tmax, tstep, tret        &
                          ,var                            &
                          ,grid                           &
#ifdef PNCREGRID
                          ,my_mpi                         &
                          ,root_pe                        &
#endif
                          )

  USE MESSY_NCREGRID_INTERFACE, ONLY: NCREGRID_SET_MESSAGEMODE &
                                    , INTERFACE_GEOHYBGRID

  IMPLICIT NONE

  ! I/O
  INTEGER, INTENT(IN)                 :: GCTRL   ! RG_SCAN, RG_PROC, RG_STOP
  INTEGER, INTENT(IN)                 :: GNML    ! NML_NEXT, NML_STAY
#ifndef PNCREGRID
  INTEGER, INTENT(OUT)                :: GSTAT   ! RGSTAT_START, RGSTAT_CONT
#else
  INTEGER, INTENT(INOUT)              :: GSTAT   ! RGSTAT_START, RGSTAT_CONT
#endif
                                                 ! RGSTAT_STOP
  CHARACTER(LEN=*) , INTENT(IN)       :: nmlfile ! namelist file

  ! OPTIONAL I/O
  INTEGER          , INTENT(IN) , OPTIONAL :: iounit   ! I/O unit for nmlfile
  LOGICAL          , INTENT(IN) , OPTIONAL :: lrgx     ! regrid along 'lon'
  LOGICAL          , INTENT(IN) , OPTIONAL :: lrgy     ! regrid along 'lat'
  LOGICAL          , INTENT(IN) , OPTIONAL :: lrgz     ! regrid along 'lev'
  INTEGER          , INTENT(IN) , OPTIONAL :: tmin     ! start time step
  INTEGER          , INTENT(IN) , OPTIONAL :: tmax     ! stop time step
  INTEGER          , INTENT(IN) , OPTIONAL :: tstep    ! time step
  INTEGER          , INTENT(IN) , OPTIONAL :: tret     ! return time step
  TYPE (ncvar), DIMENSION(:), POINTER, OPTIONAL :: var ! list of variables
#ifndef PNCREGRID
  TYPE (geohybgrid), INTENT(OUT), OPTIONAL :: grid     ! output grid info
#else
  TYPE (geohybgrid), INTENT(INOUT), OPTIONAL :: grid   ! output grid info
  TYPE (t_mpi_def) , INTENT(IN),  OPTIONAL :: my_mpi   ! MPI definition
  INTEGER          , INTENT(OUT), OPTIONAL :: root_pe  ! PE delivering data
#endif

  ! LOCAL (for DEFAULT)
  CHARACTER(LEN=*), PARAMETER :: substr = 'REGRID_CONTROL'
  INTEGER :: iou   ! I/O unit for namelist file
  LOGICAL :: lx    ! regrid along x
  LOGICAL :: ly    ! regrid along y
  LOGICAL :: lz    ! regrid along z
  ! LOCAL
  LOGICAL :: lex   ! file exists ?
  LOGICAL :: lopn  ! file opened ?
  INTEGER :: iout  ! I/O unit
#ifndef PNCREGRID
  INTEGER :: i,j   ! counter
#endif
  INTEGER :: nvars ! number of variables
  INTEGER :: mvars ! 'effective' number of variables
  INTEGER :: status
  TYPE (ncvar), DIMENSION(:), POINTER  :: xivar ! list of variables (input)
  TYPE (ncvar), DIMENSION(:), POINTER  :: xovar ! list of variables (output)
  TYPE (ncvar)                         :: tvar  ! temporal variable
  TYPE (ncvar)                         :: qvari  ! temporal variable
  TYPE (ncvar)                         :: svari  ! temporal variable (sorted)
  TYPE (ncvar)                         :: pvari  ! temporal variable (packed)
  TYPE (ncvar)                         :: qvaro  ! temporal variable
  TYPE (ncvar), DIMENSION(:), POINTER  :: svaro  ! temporal variable (sorted)
  TYPE (ncvar), DIMENSION(:), POINTER  :: pvaro  ! temporal variable (packed)

  ! NAMELIST
  INTEGER, SAVE :: fstat ! status of namelist file
  TYPE (geohybgrid)                   :: gi  ! input grid
  TYPE (geohybgrid), SAVE             :: gi0 ! input grid (raw)
  TYPE (geohybgrid)                   :: go  ! output grid
  TYPE (geohybgrid)                   :: gg  ! grdfile grid
  TYPE (geohybgrid), SAVE             :: gg0 ! grdfile grid
  CHARACTER (LEN=100*GRD_MAXSTRLEN), SAVE :: varstr  ! variable string
  CHARACTER (LEN=GRD_MAXSTRLEN), SAVE :: outfile ! netCDF output-file
  INTEGER, SAVE                 :: i_t(4) ! input time control
  INTEGER, SAVE                 :: g_t(3) ! grid time control
  INTEGER, SAVE                 :: o_t(3) ! output time control
  LOGICAL, SAVE                 :: lp     ! pressure or sigma
  LOGICAL, SAVE                 :: lint   ! input time ?
  TYPE (ncatt), SAVE            :: rggatt ! RG global attribute
  TYPE (narray), SAVE           :: nmlstr ! namelist as string
  REAL(DP), DIMENSION(2,4,2), SAVE :: ranges ! range of fields
  ! TIME CONTROL
  INTEGER                       :: ttmin, ttmax, ttstep, ttret
  INTEGER                       :: tg     ! grid time step
  INTEGER                       :: to     ! output time step
#ifndef PNCREGRID
  INTEGER                       :: tc     ! time step counter
#endif
  ! GEOHYBRID-GRIDS
  TYPE (geohybgrid)             :: gis, gix ! sorted and index input-grid
  LOGICAL                       :: ok       ! result OK ?
  TYPE (geohybgrid)             :: ggs, ggx ! sorted and index grid
  ! AXES FOR REGRIDDING
  TYPE (axis),   DIMENSION(:), POINTER :: sax, dax ! source and dest. axes
  ! N-ARRAYS FOR REGRIDDING
  TYPE (narray), DIMENSION(:), POINTER :: nai, nao ! input/output for regridder
  ! REGRIDDING STATISTICS
  REAL (DP),     DIMENSION(:), POINTER :: sovl, dovl
  INTEGER,       DIMENSION(:), POINTER :: rcnt
  ! FORCE UNLIMITED ID
  INTEGER                       :: setuid

  ! VARIABLES
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: ivar(:)   ! input variable names
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: ovar(:)   ! output variable names
  REAL(dp)                    , POINTER :: scl(:)    ! scaling
  INTEGER                     , POINTER :: RGT(:)    ! regridding type
  CHARACTER (LEN=3)           , POINTER :: RGTstr(:) ! regridding type string
  INTEGER, DIMENSION(:)       , POINTER :: dims      ! order of grid dimensions
  INTEGER                               :: axes(3)   ! grid axes of variable
#ifdef PNCREGRID
  INTEGER, SAVE       :: GSTAT_SAVE

  ! If running in MPI environment, the variables can be handled in parallel.
  ! Parallel processing is enabled by the presence of my_mpi.

  IF(PRESENT(my_mpi))  THEN
     me    = my_mpi%rank
     MY_RANK = me
     nproc = my_mpi%nproc
     comm  = my_mpi%comm
  ELSE
     me    = 0
     nproc = 1
     comm  = -1
     work_loop = 0
  END IF

  ! For every namelist, REGRID_CONTROL has to be called n cycles.
  ! During the first iteration of a cycle, all work is done in parallel.
  ! During the next iterations, one PE is returned as root_pe to enable
  ! the calling routine to process the data computed on root_pe.
 
  IF(work_loop == 0)   THEN    ! all worker PEs compute in parallel
     CALL REGRID_CONTROL_INIT

     IF(i_am_worker)   THEN
        CALL REGRID_CONTROL_WORK
     END IF

     GSTAT_SAVE = GSTAT
  END IF

  IF(PRESENT(root_pe)) THEN
     root_pe = pe_list(work_loop)
  END IF

  IF(work_loop == nproc_work-1)  THEN
     work_loop = 0
     GSTAT = GSTAT_SAVE
  ELSE
     work_loop = work_loop+1
     IF(GSTAT_SAVE == RGSTAT_STOP)   THEN   ! no cycle, if all work is done
        GSTAT = RGSTAT_STOP
     ELSE
        GSTAT = RGSTAT_CYCLE
     END IF
  END IF

  RETURN

 CONTAINS

 SUBROUTINE REGRID_CONTROL_INIT

   IMPLICIT NONE

   INTEGER :: i
#endif

  ! INIT
  CALL NCREGRID_SET_MESSAGEMODE

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
  IF (PRESENT(lrgx)) THEN
     lx = lrgx
  ELSE
     lx = .TRUE.               ! DEFAULT
  END IF
  !
  IF (PRESENT(lrgy)) THEN
     ly = lrgy
  ELSE
     ly = .TRUE.               ! DEFAULT
  END IF
  !
  IF (PRESENT(lrgz)) THEN
     lz = lrgz
  ELSE
     lz = .TRUE.               ! DEFAULT
  END IF
  !
  ! INITIALIZE OUTPUT PARAMETER
  IF (PRESENT(var)) THEN
     IF (ASSOCIATED(var)) THEN
        DO i=1, SIZE(var)
           CALL INIT_NCVAR(var(i))
        END DO
        DEALLOCATE(var, STAT=status)
        CALL ERRMSG(substr,status,1)
     END IF
     NULLIFY(var)
  END IF
  IF (PRESENT(grid)) CALL INIT_GEOHYBGRID(grid)

  ! INIT
  NULLIFY(dims)
  NULLIFY(ivar)
  NULLIFY(ovar)
  NULLIFY(scl)
  NULLIFY(RGT)
  NULLIFY(RGTstr)
  NULLIFY(sax)
  NULLIFY(dax)
  NULLIFY(xivar)
  NULLIFY(xovar)
  NULLIFY(svaro)
  NULLIFY(pvaro)
  NULLIFY(nai)
  NULLIFY(nao)
  NULLIFY(sovl)
  NULLIFY(dovl)
  NULLIFY(rcnt)

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
     CALL RGMSG(substr, RGMLVL, &
         '-------------------------------------------------------')
     CALL RGMSG(substr, RGMLVL, &
          'START REGRIDDING PROCEDURE (NCREGRID VERSION '//&
          &TRIM(NCREGRIDVERS)//')' )
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
     CALL RGMSG(substr, RGMLI, &
          'READING NAMELIST FROM FILE '''//TRIM(nmlfile)//''' ...')
     CALL READ_NAMELIST(iou, fstat               &
                        ,gi0, gg0, outfile       &
                        ,varstr                  &
                        ,i_t, g_t, o_t, lp, lint, nmlstr, ranges)
     ! CHECK END OF FILE
     IF (fstat == 0) THEN    ! GO ON
        CALL RGMSG(substr, RGMLIC, '... OK !')
#ifdef PNCREGRID
        ! NOTE: THE FOLLOWING SUBROUTINE CALLS ARE REQUIRED TO SET
        !       i_am_worker CORRECTLY
        ! GET VARIABLE INFORMATION
        CALL PARSE_VARSTR(varstr, ivar, ovar, scl     &
             ,RGT, RGTstr,TRIM(gi0%file))
        CALL DISTRIBUTE_VARS_ON_PES(ivar, ovar, scl, RGT, RGTstr)
#endif
        GSTAT = RGSTAT_CONT
     ELSE                    ! ERROR OR END OF FILE
#ifdef PNCREGRID
        i_am_worker = .FALSE.
#endif
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
  END IF

  IF (GCTRL == RG_STOP) GSTAT = RGSTAT_STOP

  ! --------------------------------------------------
  IF (GSTAT == RGSTAT_STOP) THEN
     CLOSE(iou)
     ! CLEAN
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
           CALL ERRMSG(substr,status,2)
        END IF
        NULLIFY(var)
     END IF
     !
     CALL INIT_GEOHYBGRID(gi0)
     CALL INIT_GEOHYBGRID(gg0)
     CALL INIT_GEOHYBGRID(go)
     CALL INIT_NARRAY(nmlstr)
     CALL INIT_NCATT(rggatt)
     !
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

#ifdef PNCREGRID
  RETURN

  END SUBROUTINE REGRID_CONTROL_INIT

  SUBROUTINE REGRID_CONTROL_WORK

    IMPLICIT NONE

    INTEGER :: i, j, tc

  IF (GSTAT == RGSTAT_STOP) RETURN
#endif

  ! --- TIME CONTROL ---------------------------------
  ! UPDATE TIME CONTROL FROM NAMELIST
  ttmin = i_t(1)
  ttstep = i_t(2)
  ttmax = i_t(3)
  ttret = i_t(4)
  !
  ! OVERWRITE TIME CONTROL FROM SUBROUTINE CALL
  IF (PRESENT(tmin)) THEN
     ttmin = tmin
  END IF
  IF (PRESENT(tmax)) THEN
     ttmax = tmax
  END IF
  IF (PRESENT(tstep)) THEN
     ttstep = tstep
  END IF
  IF (PRESENT(tret)) THEN
     ttret = tret
  END IF
  !
  ! CHECK TIME CONTROL SETTINGS
  ! CHECK INFILE TIME STEPPING / INTERFACE MODE TIME STEPPING
  IF (ttmin <= 0) THEN
     CALL RGMSG(substr, RGMLE,  'i_t(1) IN NAMELIST OR tmin AT ', .false.)
     CALL RGMSG(substr, RGMLEC, 'SUBROUTINE CALL MUST BE >0 !')
  END IF
  IF (ttstep <= 0) THEN
     CALL RGMSG(substr, RGMLE,  'i_t(2) IN NAMELIST OR tstep AT ', .false.)
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
     CALL RGMSG(substr, RGMLEC, 'SUBROUTINE CALL MUST BE >= i_t(1) OR tmin!')
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
  ! NOTE: IN PRINCIPLE IT IS POSSIBLE TO SPECIFY A GRDFILE (IN THE NAMELIST)
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

  ! GET VARIABLE INFORMATION
  CALL PARSE_VARSTR(varstr, ivar, ovar, scl     &
                    ,RGT, RGTstr,TRIM(gi0%file))
#ifdef PNCREGRID
  CALL RGMSG(substr, RGMLI, &
       'DISTRIBUTING VARIABLES ON PEs ...')
  CALL DISTRIBUTE_VARS_ON_PES(ivar, ovar, scl, RGT, RGTstr)
  DO i=1, SIZE(ivar)
     IF (i_am_worker) &
          CALL RGMSG(substr, RGMLIC &
          , '... ... '//TRIM(ivar(i))//' on PE=', me, '')
  END DO
  !
  CALL EXPAND_FILENAME(outfile, SIZE(ivar), me)
  IF (i_am_worker) & 
       CALL RGMSG(substr, RGMLIC, '... OUTPUT FILE '//TRIM(outfile))
  CALL RGMSG(substr, RGMLIC, '... OK !')
#endif
  nvars = SIZE(ivar)
  ALLOCATE(xivar(nvars), STAT=status)
  CALL ERRMSG(substr,status,3)

  ! START TIME LOOP
  tg = g_t(1) - g_t(2)              ! INITIALIZE time step in grid file
  to = o_t(1) - o_t(2)              ! INITIALIZE time step in output file
  DO tc = ttmin,ttmax, ttstep
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
     gi%t = tc
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
     CALL CHECK_GEOHYBGRID(gi, ranges(1,:,:))
     CALL RGMSG(substr, RGMLIC, ' ... O.K.')

     ! GET DESTINATION GRID
     CALL RGMSG(substr, RGMLVL, 'IMPORTING OUTPUT GRID ...')
     IF (TRIM(gg%file) == '') THEN
        CALL RGMSG(substr, RGMLVM, '... CHECKING FOR INTERFACE ...')
        CALL INTERFACE_GEOHYBGRID(gg, ok)
        IF (.NOT.ok) THEN
           CALL RGMSG(substr, RGMLE, &
                'EITHER GRDFILE HAS TO BE SPECIFIED IN NAMELIST,' &
                , .false.)
           CALL RGMSG(substr, RGMLEC, &
                'OR INTERFACE HAS TO BE COMPILED INTO NCREGRID !')
        END IF
        IF (lint) THEN  ! ADJUST OUTPUT TIME STEPPING FOR INTERFACE MODE
           to = tc
           go%t = gi%t
        ELSE
           to = tg
           go%t = gg%t
        END IF
     ELSE
        CALL RGMSG(substr, RGMLVL, &
             '... IMPORTING FROM '''//TRIM(gg%file)//'''')
        CALL IMPORT_GEOHYBGRID(gg)
     END IF
     !
     CALL RGMSG(substr, RGMLIC, ' ... SWITCHING ...')
     CALL SWITCH_GEOHYBGRID(gg, lx, ly, lz)
     CALL RGMSG(substr, RGMLIC, ' ... CHECKING ...')
     CALL CHECK_GEOHYBGRID(gg, ranges(2,:,:))
     CALL RGMSG(substr, RGMLIC, ' ... O.K.')

     ! TIME BALANCING
     CALL RGMSG(substr, RGMLVM, '>>> BALANCING TIME AXIS ...')
     CALL BALANCE_GEOHYBGRID_TIME(gi, gg, lint)
     CALL RGMSG(substr, RGMLVM, '<<< ... END BALANCING TIME AXIS !')

     IF (gi%col%dat%type == VTYPE_UNDEF) THEN
        ! SURFACE PRESSURE
        CALL RGMSG(substr, RGMLVM, '>>> BALANCING SURFACE PRESSURE ...')
        CALL BALANCE_GEOHYBGRID_PS(gi, gg, ranges)
        CALL RGMSG(substr, RGMLVM, '<<< ... END BALANCING SURFACE PRESSURE !')

        CALL RGMSG(substr, RGMLVM, '>>> SORTING INPUT GRID ...')
        CALL SORT_GEOHYBGRID(gi, gis, gix)
        CALL RGMSG(substr, RGMLVM, '<<< ... END SORTING INPUT GRID !')

        CALL RGMSG(substr, RGMLVM, '>>> SORTING OUTPUT GRID ...')
        CALL SORT_GEOHYBGRID(gg, ggs, ggx)
        CALL RGMSG(substr, RGMLVM, '<<< ... END SORTING OUTPUT GRID !')

        CALL RGMSG(substr, RGMLVM, '>>> COMPLETING INPUT GRID ...')
        CALL COMPLETE_GEOHYBGRID(gis, ranges(1,:,:), gix)
        CALL RGMSG(substr, RGMLVM, '<<< ... END COMPLETING INPUT GRID !')

        CALL RGMSG(substr, RGMLVM, '>>> COMPLETING OUTPUT GRID ...')
        CALL COMPLETE_GEOHYBGRID(ggs, ranges(2,:,:), ggx)
        CALL RGMSG(substr, RGMLVM, '<<< ... END COMPLETING OUTPUT GRID !')

        ! BALANCE GRID
        CALL RGMSG(substr, RGMLVM, '>>> BALANCING INPUT/OUTPUT GRID ...')
        CALL BALANCE_GEOHYBGRID(gis, ggs)
        CALL BALANCE_GEOHYBGRID(gix, ggx)
        ! SAVE OUTPUT TIME STEP IN GRID-STRUCTURE
        ggs%t = to  ! OUTPUT TIME STEP
        ggx%t = to  ! OUTPUT TIME STEP
        CALL RGMSG(substr, RGMLVM, '<<< END BALANCING INPUT/OUTPUT GRID !')

        ! CONVERT HYBRID GRIDS TO AXES LIST PAIRS
        CALL RGMSG(substr, RGMLVM, '>>> CONSTRUCTING REGRIDDING AXES ...')
        CALL GEOHYBGRID_AXES(gis, sax, ggs, dax, lp)
        CALL RGMSG(substr, RGMLVM, '<<< END CONSTRUCTING REGRIDDING AXES !')
     END IF

     ! GET INPUT VAR
     CALL RGMSG(substr, RGMLVL, &
          'IMPORTING ',nvars,' VARIABLE(S) FROM FILE '''//TRIM(gi%file)//'''')
     CALL RGMSG(substr, RGMLVL, &
          '('''//TRIM(gi%timem%name)//''' STEP: ',tc,')')
     mvars = 0
     DO i=1, nvars  ! LOOP OVER VARIABLE NAMES
        CALL RGMSG(substr, RGMLVL,' ... '''//TRIM(ivar(i))//''' ...')
        IF (ASSOCIATED(gi%timem%dim)) THEN
           setuid = gi%timem%dim(1)%id
        ELSE
           setuid = NULL_DIMID
        END IF
        CALL IMPORT_NCVAR(tvar, tc, varname=TRIM(ivar(i)) &
                          ,file=TRIM(gi%file)             &
                          ,setuid=setuid)
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
           ! RENAME
           IF (TRIM(ivar(i)) /= TRIM(ovar(i))) THEN
              CALL RGMSG(substr, RGMLVLC, &
                   '  ... RENAMING TO '''//TRIM(ovar(i))//'''... ')
              CALL RENAME_NCVAR(tvar, TRIM(ovar(i)))
              CALL ADD_NCATT(tvar, 'RG_ORIGINAL_NAME', vs=TRIM(ivar(i)))
           END IF
           ! SCALE
           IF (ABS(scl(i) - 1.0_dp) >= TINY(1.0_dp)) THEN
              CALL  RGMSG(substr, RGMLVLC,'  ... SCALING BY ',scl(i),' ...')
              CALL SCALE_NARRAY(tvar%dat, scl(i))
              CALL ADD_NCATT(tvar, 'RG_SCALED_BY',      &
                   vd=(/ scl(i) /))
           END IF
           ! ADJUST TIME AXIS
           ! if input_time (= lint) == .false.
           ! variable must get correct time-dimension-ID
           IF (.NOT.lint) THEN
              IF (.NOT.(setuid == NULL_DIMID)) THEN
                 DO j=1, tvar%ndims
!                    IF (tvar%dim(j)%id == setuid) THEN
                    IF (tvar%dim(j)%fuid) THEN
                       EXIT
                    END IF
                 END DO
                 CALL RGMSG(substr, RGMLVMC,'  ... ADJUSTING TIME-DIM.: ')
                 CALL RGMSG(substr, RGMLVMC,'      '//TRIM(tvar%dim(j)%name)&
                      //' (ID = ',j,') -> ')
                 CALL RGMSG(substr, RGMLVMC,'      '&
                      &//TRIM(gg%timem%dim(1)%name)&
                      //' (ID = ',setuid,') ...')
                 CALL COPY_NCDIM(tvar%dim(j), gg%timem%dim(1))
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
        IF (tc == ttret) THEN      ! SELECTED TIME STEP
           CALL RGMSG(substr, RGMLVM, 'TIME STEP ...',ttret,' ')
           IF (PRESENT(var)) THEN
              CALL RGMSG(substr, RGMLVM, '... RETURNING RAW DATA')
              ALLOCATE(var(mvars), STAT=status)
              CALL ERRMSG(substr,status,5)
              DO i=1, mvars
                 CALL INIT_NCVAR(var(i))
                 CALL COPY_NCVAR(var(i),xivar(i))
              END DO
           END IF
           IF (PRESENT(grid)) THEN
              CALL RGMSG(substr, RGMLVM, '... RETURNING INPUT GRID')
              CALL COPY_GEOHYBGRID(grid, gi)
           END IF
        END IF
        CALL RGMSG(substr, RGMLVL, '... DONE (SCAN MODE) !')
     ! ---------------------------------------------------
     CASE(RG_PROC)                ! REGRID
     ! ---------------------------------------------------
        CALL RGMSG(substr, RGMLVL, 'REGRIDDING MODE ...')
        ! ALLOCATE SPACE
        ALLOCATE(xovar(mvars), svaro(mvars), pvaro(mvars), STAT=status)
        CALL ERRMSG(substr,status,6)
        ALLOCATE(nai(mvars), STAT=status)
        CALL ERRMSG(substr,status,7)
        DO i=1, mvars             ! LOOP OVER ALL VALID VARIABLES
           CALL RGMSG(substr, RGMLVL, &
                'VARIABLE '''//TRIM(xivar(i)%name)//''' ...')
           ! CONVERT IDX fields to FRACTION
           IF ((RGT(i) == RG_IDX).OR.(RGT(i) == RG_IXF)) THEN
              CALL RGMSG(substr, RGMLVMC, ' ... INDEX FRACTIONS ...')
              CALL IDX2FRAC_NCVAR(xivar(i), qvari)
              CALL RGMSG(substr, RGMLVMC, ' ... ->'''//TRIM(qvari%name)//'''')
           ELSE
              CALL COPY_NCVAR(qvari, xivar(i))
           END IF
           ! GET ORDER INFORMATION
           CALL RGMSG(substr, RGMLIC, ' ... CHECKING ...')
           CALL CHECK_NCVAR_ON_GEOHYBGRID(qvari, gis, dims, axes, ok)
           ! BALANCE OUTPUT VARIABLE
           CALL RGMSG(substr, RGMLIC, ' ... BALANCING ...')
           CALL BALANCE_GEOHYBGRID_NCVAR(qvari, axes, ggs, xovar(i))
           ! SORT VARIABLE ACCORDING TO GRID
           CALL RGMSG(substr, RGMLIC, ' ... SORTING ...')
           CALL SORT_GEOHYBGRID_NCVAR(qvari, gix, axes, svari)
           CALL BALANCE_GEOHYBGRID_NCVAR(svari, axes, ggs, svaro(i))
           ! PACK VARIABLE
           CALL RGMSG(substr, RGMLIC, ' ... PACKING ...')
           CALL PACK_GEOHYBGRID_NCVAR(svari, dims, axes, pvari)
           CALL BALANCE_GEOHYBGRID_NCVAR(pvari, axes, ggs, pvaro(i))
           CALL INIT_NARRAY(pvaro(i)%dat)
           ! FILL N-ARRAY
           CALL COPY_NARRAY(nai(i), pvari%dat)
           ! CLEAN
           CALL INIT_NCVAR(qvari)
           CALL INIT_NCVAR(svari)
           CALL INIT_NCVAR(pvari)
           DEALLOCATE(dims, STAT=status)
           CALL ERRMSG(substr,status,8)
           NULLIFY(dims)
           axes(:) = 0
           CALL RGMSG(substr, RGMLIC, ' ... DONE (VARIABLE) !')
        END DO
        !
        ! REGRID
        IF (mvars > 0) THEN
           CALL NREGRID(nai, sax, dax, nao, RGT(1:mvars), sovl, dovl, rcnt)
           ! OUTPUT STATISTICS
           IF (MSGMODE > MSGMODE_S) THEN
              CALL NREGRID_STAT(sax, dax, sovl, dovl, nai, nao, rcnt)
           END IF
           !
           ! CLEAN
           DEALLOCATE(sovl, STAT=status)
           CALL ERRMSG(substr,status,9)
           NULLIFY(sovl)
           DEALLOCATE(dovl, STAT=status)
           CALL ERRMSG(substr,status,10)
           NULLIFY(dovl)
           DEALLOCATE(rcnt, STAT=status)
           CALL ERRMSG(substr,status,11)
           NULLIFY(rcnt)
        END IF
        !
        !
        DO i=1, mvars
           CALL RGMSG(substr, RGMLVL, &
                'VARIABLE '''//TRIM(xovar(i)%name)//''' ...')
           ! COPY N-ARRAY
           CALL COPY_NARRAY(pvaro(i)%dat, nao(i))   ! 'DATA'
           CALL CHECK_NCVAR_ON_GEOHYBGRID(svaro(i), ggs, dims, axes, ok)
           ! UNPACK VARIABLE
           CALL RGMSG(substr, RGMLIC, ' ... UN-PACKING ...')
           CALL PACK_GEOHYBGRID_NCVAR(pvaro(i), dims, axes, svaro(i), .true.)
           !  UN-SORT
           CALL RGMSG(substr, RGMLIC, ' ... UN-SORTING ...')
           CALL SORT_GEOHYBGRID_NCVAR(svaro(i), ggx, axes, qvaro, .true.)
           !  UN-IDX
           CALL INIT_NCVAR(xovar(i))
           IF (RGT(i) == RG_IDX) THEN
              CALL RGMSG(substr, RGMLVMC, ' ... MAXIMUM INDEX FRACTION ...')
              CALL MAXFRAC2IDX_NCVAR(qvaro, xovar(i))
              xovar(i)%name = TRIM(xivar(i)%name)
              CALL RGMSG(substr, RGMLVMC, &
                   ' ... ->'''//TRIM(xovar(i)%name)//'''')
           ELSE
              CALL COPY_NCVAR(xovar(i), qvaro)
           END IF
           ! CLEAN
           CALL INIT_NCVAR(pvaro(i))
           CALL INIT_NCVAR(svaro(i))
           CALL INIT_NCVAR(qvaro)
           DEALLOCATE(dims, STAT=status)
           CALL ERRMSG(substr,status,12)
           NULLIFY(dims)
           axes(:) = 0
           CALL RGMSG(substr, RGMLIC, ' ... DONE (VARIABLE) !')
        END DO
        ! CLEAN
        DO i=1, mvars
           CALL INIT_NARRAY(nai(i))
           CALL INIT_NARRAY(nao(i))
        END DO
        IF (mvars > 0) THEN
           DEALLOCATE(nai, STAT=status)
           CALL ERRMSG(substr,status,13)
           NULLIFY(nai)
           DEALLOCATE(nao, STAT=status)
           CALL ERRMSG(substr,status,14)
           NULLIFY(nao)
           DEALLOCATE(svaro, pvaro, STAT=status)
           CALL ERRMSG(substr,status,15)
           NULLIFY(svaro)
           NULLIFY(pvaro)
        END IF
        !
        ! UN-SORT COMPLETED AND BALANCED GRID
        CALL SORT_GEOHYBGRID(ggs, go, ggx, .true.)
        !
        ! OUTPUT GRID AND VARIABLES
        IF (TRIM(outfile) /= '') THEN
           go%file = TRIM(outfile)    ! name of output-file
           go%t = to                  ! output time step
           CALL EXPORT_GEOHYBGRID(go)
           ! EXPORT GLOBL ATTRIBUTE, ONLY AT FIRST TIME STEP
           IF (tc == ttmin) THEN
              ! GET OLD ATTRIBUTE
              CALL IMPORT_NCATT(rggatt, varid=NF90_GLOBAL   &
                                ,attname = 'RG_HISTORY'     &
                                ,file = TRIM(go%file), lnostop=.true.)
              ! CHECK IF EMPTY OR CHAR-ATT
              IF ((rggatt%dat%type == VTYPE_UNDEF).OR.  &
                  (rggatt%dat%type == VTYPE_CHAR)) THEN
                 ! APPEND NEW DATA
                 CALL CAT_NARRAY(rggatt%dat, nmlstr)
                 rggatt%xtype = NF90_CHAR
                 rggatt%len   = rggatt%dat%dim(1)
              ELSE  ! EXISTS WITH NON-CHAR TYPE
                 CALL RGMSG(substr, RGMLW, &
                      'ATTRIBUTE '''//TRIM(rggatt%name)//'''')
                 CALL RGMSG(substr, RGMLWC, &
                      'EXISTS ALREADY, BUT IS NOT OF TYPE CHAR !')
              END IF
              ! WRITE ATTRIBUTE
              CALL EXPORT_NCATT(rggatt, file=TRIM(go%file)  &
                                ,clobber=.true.)
           END IF
           ! EXPORT VARIABLES
           DO i=1, mvars
              xovar(i)%ustep = go%t
              CALL EXPORT_NCVAR(xovar(i), file=TRIM(go%file))
           END DO
        END IF
        !
        ! RETURN VALUES TO SUBROUTINE CALL
        IF (tc == ttret) THEN      ! SELECTED TIME STEP
           CALL RGMSG(substr, RGMLVM, &
                'RETURN STEP REACHED: TIME STEP ',ttret,'')
           IF (PRESENT(var)) THEN
              CALL RGMSG(substr, RGMLVM, '... RETURNING REGRIDDED DATA')
              ALLOCATE(var(mvars), STAT=status)
              CALL ERRMSG(substr,status,16)
              DO i=1, mvars
                 xovar(i)%ustep = go%t
                 CALL COPY_NCVAR(var(i), xovar(i))
              END DO
           END IF
           IF (PRESENT(grid)) THEN
              CALL RGMSG(substr, RGMLVM, '... RETURNING OUTPUT GRID')
              CALL COPY_GEOHYBGRID(grid, go)
           END IF
        END IF
        !
        ! CLEAN
        DO i=1, mvars
           CALL INIT_NCVAR(xovar(i))
        END DO
        DEALLOCATE(xovar, STAT=status)
        CALL ERRMSG(substr,status,17)
        NULLIFY(xovar)
        !
        CALL RGMSG(substr, RGMLVL, '... DONE (REGRIDDING MODE) !')
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
     IF (ASSOCIATED(sax)) THEN
        DO i=1, SIZE(sax)
           CALL INIT_AXIS(sax(i))
        END DO
        DEALLOCATE(sax, STAT=status)
        CALL ERRMSG(substr,status,18)
        NULLIFY(sax)
     END IF
     !
     IF (ASSOCIATED(dax)) THEN
        DO i=1, SIZE(dax)
           CALL INIT_AXIS(dax(i))
        END DO
        DEALLOCATE(dax, STAT=status)
        CALL ERRMSG(substr,status,19)
        NULLIFY(dax)
     END IF
     !
     CALL INIT_GEOHYBGRID(gis)
     CALL INIT_GEOHYBGRID(gix)
     CALL INIT_GEOHYBGRID(gg)
     CALL INIT_GEOHYBGRID(ggx)
     CALL INIT_GEOHYBGRID(ggs)
     CALL INIT_GEOHYBGRID(gi)
     CALL INIT_GEOHYBGRID(go)

  END DO  ! TIME LOOP

  ! CLEAN
  IF (ASSOCIATED(ivar)) DEALLOCATE(ivar)
  IF (ASSOCIATED(ovar)) DEALLOCATE(ovar)
  IF (ASSOCIATED(scl))  DEALLOCATE(scl)
  IF (ASSOCIATED(RGT))  DEALLOCATE(RGT)
  IF (ASSOCIATED(RGTstr)) DEALLOCATE(RGTstr)
  NULLIFY(ivar, ovar, scl, RGT, RGTstr)

  IF (ASSOCIATED(xivar)) THEN
     ! CALL INIT_NCVAR(xivar(i)) HAS ALREADY BEEN DONE INSIDE
     ! (AT THE END OF) THE TIME LOOP
     DEALLOCATE(xivar, STAT=status)
     CALL ERRMSG(substr,status,20)
     NULLIFY(xivar)
  END IF

#ifdef PNCREGRID
  END SUBROUTINE REGRID_CONTROL_WORK
#endif

END SUBROUTINE REGRID_CONTROL
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE READ_NAMELIST(iou, status                    &
                         ,gi, gg, outfile               &
                         ,var, i_t, g_t, o_t, lp, lint  &
                         ,nmlstr, ranges)

  IMPLICIT NONE

  ! I/O
  INTEGER                      , INTENT(IN)  :: iou    ! logical I/O unit
  INTEGER                      , INTENT(OUT) :: status ! file status
  TYPE (geohybgrid)            , INTENT(OUT) :: gi       ! input grid
  TYPE (geohybgrid)            , INTENT(OUT) :: gg       ! output grid
  CHARACTER (LEN=GRD_MAXSTRLEN), INTENT(OUT) :: outfile  ! output filename
  CHARACTER (LEN=100*GRD_MAXSTRLEN), INTENT(OUT) :: var  ! variable string
  INTEGER                      , INTENT(OUT) :: i_t(4)   ! input time control
  INTEGER                      , INTENT(OUT) :: g_t(3)   ! grid time control
  INTEGER                      , INTENT(OUT) :: o_t(3)   ! output time control
  LOGICAL                      , INTENT(OUT) :: lp       ! pressure or sigma
  LOGICAL                      , INTENT(OUT) :: lint     ! input time ?
  TYPE (narray)                , INTENT(OUT) :: nmlstr   ! namelist as string
  REAL(DP), DIMENSION(2,4,2)   , INTENT(OUT) :: ranges   ! range of fields

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'READ_NAMELIST'
  CHARACTER(GRD_MAXSTRLEN) :: infile            ! input netCDF filename
! CHARACTER(GRD_MAXSTRLEN) :: outfile           ! output netCDF filename
  CHARACTER(GRD_MAXSTRLEN) :: grdfile           ! grid file (netCDF)
  CHARACTER(GRD_MAXSTRLEN) :: i_latm, i_lonm    ! input latitude/longitude
  CHARACTER(GRD_MAXSTRLEN) :: i_lati, i_loni    ! input latitude/longitude
  CHARACTER(GRD_MAXSTRLEN) :: i_col             ! input columns
  CHARACTER(GRD_MAXSTRLEN) :: g_latm, g_lonm    ! output latitude/longitude
  CHARACTER(GRD_MAXSTRLEN) :: g_lati, g_loni    ! output latitude/longitude
  CHARACTER(GRD_MAXSTRLEN) :: g_col             ! output columns
  CHARACTER(GRD_MAXSTRLEN) :: i_hyai, i_hybi    ! input hybrid coeff.
  CHARACTER(GRD_MAXSTRLEN) :: i_hyam, i_hybm    ! input hybrid coeff. (mid)
  CHARACTER(GRD_MAXSTRLEN) :: g_hyai, g_hybi    ! output hybrid coeff.
  CHARACTER(GRD_MAXSTRLEN) :: g_hyam, g_hybm    ! output hybrid coeff. (mid)
  CHARACTER(GRD_MAXSTRLEN) :: i_timei, i_timem  ! input time
  CHARACTER(GRD_MAXSTRLEN) :: g_timei, g_timem  ! output time
  CHARACTER(GRD_MAXSTRLEN) :: i_ps, g_ps        ! input/output surface press.
  CHARACTER(GRD_MAXSTRLEN) :: i_p0, g_p0        ! input/output ref. press.
  LOGICAL                  :: pressure          ! pressure or sigma
  LOGICAL                  :: input_time        ! input time information
  !
  REAL, DIMENSION(2)       :: i_latr, i_lonr
  REAL, DIMENSION(2)       :: i_hyar, i_hybr
  REAL, DIMENSION(2)       :: g_latr, g_lonr
  REAL, DIMENSION(2)       :: g_hyar, g_hybr


  ! HELPERS
  CHARACTER(GRD_MAXSTRLEN) :: tstr
  CHARACTER(LEN=8)         :: date
  CHARACTER(LEN=10)        :: time
  CHARACTER(LEN=5)         :: zone

  NAMELIST /regrid/ infile, outfile, grdfile                           &
       ,i_latm, i_lati, i_lonm, i_loni, i_hyam, i_hyai, i_hybm, i_hybi &
       ,i_timei, i_timem, i_ps, i_p0                                   &
       ,i_col                                                          &
       ,g_latm, g_lati, g_lonm, g_loni, g_hyam, g_hyai, g_hybm, g_hybi &
       ,g_timei, g_timem, g_ps, g_p0                                   &
       ,g_col                                                          &
       ,var, i_t, g_t, o_t, pressure, input_time                       &
       ,i_latr, i_lonr, i_hyar, i_hybr                                 &
       ,g_latr, g_lonr, g_hyar, g_hybr

  ! INIT (DEFAULT VALUES)
  infile = ''
  outfile = ''
  grdfile = ''
  i_latm = ''
  i_lati = ''
  i_lonm = ''
  i_loni = ''
  i_col = ''
  g_col = ''
  i_hyam = ''
  i_hyai = ''
  i_hybm = ''
  i_hybi = ''
  i_timei = ''
  i_timem = ''
  i_ps = ''
  i_p0 = ''
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
  var = ''
  i_t(:) = (/ 1,1,0,0 /)
  g_t(:) = (/ 1,1,0 /)
  o_t(:) = (/ 1,1,0 /)
  pressure   = .FALSE.  ! DEFAULT IS: SIGMA LEVELS
  input_time = .TRUE.   ! DEFAULT IS: USE INPUT TIME INFO FOR OUTPUT

  ranges(:,:,:) = RGEMPTY
  i_latr(:) = RGEMPTY
  i_lonr(:) = RGEMPTY
  i_hyar(:) = RGEMPTY
  i_hybr(:) = RGEMPTY
  g_latr(:) = RGEMPTY
  g_lonr(:) = RGEMPTY
  g_hyar(:) = RGEMPTY
  g_hybr(:) = RGEMPTY

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
  CALL ADDLINE_CNARRAY(nmlstr, 'NCREGRID VERSION '//TRIM(NCREGRIDVERS))
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
  CALL ADDLINE_CNARRAY(nmlstr, 'i_col   = '''//TRIM(i_col)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'g_col   = '''//TRIM(g_col)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_hyam  = '''//TRIM(i_hyam)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_hyai  = '''//TRIM(i_hyai)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_hybm  = '''//TRIM(i_hybm)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_hybi  = '''//TRIM(i_hybi)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_timei = '''//TRIM(i_timei)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_timem = '''//TRIM(i_timem)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_ps    = '''//TRIM(i_ps)//''',')
  CALL ADDLINE_CNARRAY(nmlstr, 'i_p0    = '''//TRIM(i_p0)//''',')
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
  WRITE(tstr,'(e10.4,",",e10.4)') i_lonr(1),i_lonr(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'i_lonr  = '//TRIM(tstr))
  WRITE(tstr,'(e10.4,",",e10.4)') i_latr(1),i_latr(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'i_latr  = '//TRIM(tstr))
  WRITE(tstr,'(e10.4,",",e10.4)') i_hyar(1),i_hyar(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'i_hyar  = '//TRIM(tstr))
  WRITE(tstr,'(e10.4,",",e10.4)') i_hybr(1),i_hybr(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'i_hybr  = '//TRIM(tstr))
  !
  WRITE(tstr,'(e10.4,",",e10.4)') g_lonr(1),g_lonr(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'g_lonr  = '//TRIM(tstr))
  WRITE(tstr,'(e10.4,",",e10.4)') g_latr(1),g_latr(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'g_latr  = '//TRIM(tstr))
  WRITE(tstr,'(e10.4,",",e10.4)') g_hyar(1),g_hyar(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'g_hyar  = '//TRIM(tstr))
  WRITE(tstr,'(e10.4,",",e10.4)') g_hybr(1),g_hybr(2)
  CALL ADDLINE_CNARRAY(nmlstr, 'g_hybr  = '//TRIM(tstr))
  !
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
  gi%col%name   = TRIM(i_col)
  gi%hyam%name  = TRIM(i_hyam)
  gi%hyai%name  = TRIM(i_hyai)
  gi%hybm%name  = TRIM(i_hybm)
  gi%hybi%name  = TRIM(i_hybi)
  gi%timei%name = TRIM(i_timei)
  gi%timem%name = TRIM(i_timem)
  gi%ps%name    = TRIM(i_ps)
  gi%p0%name    = TRIM(i_p0)

  gg%file       = TRIM(grdfile)
  gg%latm%name  = TRIM(g_latm)
  gg%lati%name  = TRIM(g_lati)
  gg%lonm%name  = TRIM(g_lonm)
  gg%loni%name  = TRIM(g_loni)
  gg%col%name   = TRIM(g_col)
  gg%hyam%name  = TRIM(g_hyam)
  gg%hyai%name  = TRIM(g_hyai)
  gg%hybm%name  = TRIM(g_hybm)
  gg%hybi%name  = TRIM(g_hybi)
  gg%timei%name = TRIM(g_timei)
  gg%timem%name = TRIM(g_timem)
  gg%ps%name    = TRIM(g_ps)
  gg%p0%name    = TRIM(g_p0)

  lp   = pressure
  lint = input_time

!  var
!  i_t
!  g_t
!  o_t
!  outfile

  ranges(1,1,:) = i_lonr(:)
  ranges(1,2,:) = i_latr(:)
  ranges(1,3,:) = i_hyar(:)
  ranges(1,4,:) = i_hybr(:)
  ranges(2,1,:) = g_lonr(:)
  ranges(2,2,:) = g_latr(:)
  ranges(2,3,:) = g_hyar(:)
  ranges(2,4,:) = g_hybr(:)

END SUBROUTINE READ_NAMELIST
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE ADDLINE_CNARRAY(na, str)

  IMPLICIT NONE

  ! I/O
  TYPE (narray),    INTENT(INOUT) :: na
  CHARACTER(LEN=*), INTENT(IN)    :: str

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'ADDLINE_CNARRAY'
  INTEGER :: n, m ! string length
  INTEGER :: i ! counter
  TYPE (narray) :: nas

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

! ------------------------------------------------------------------
SUBROUTINE PARSE_VARSTR(var, ivar, ovar, scl, RGT, RGTstr, file)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)           :: var       ! variable string
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER  :: ivar(:)   ! input variable names
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER  :: ovar(:)   ! output variabel names
  REAL(dp)                    , POINTER  :: scl(:)    ! scaling factor
  INTEGER                     , POINTER  :: RGT(:)    ! regridding type
  CHARACTER(LEN=3)            , POINTER  :: RGTstr(:) ! regridding type string
  CHARACTER(LEN=*)            , OPTIONAL, INTENT(IN) :: file      ! filename

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'PARSE_VARSTR'
  TYPE (ncvar), DIMENSION(:), POINTER   :: kvar   ! variables from file
  INTEGER                               :: nvar   ! number of variables
  INTEGER                               :: status
  INTEGER                               :: i,j
  INTEGER                               :: idx1, idx2, idx3
  INTEGER                               :: idxl, idxr
  CHARACTER(GRD_MAXSTRLEN), ALLOCATABLE :: lvar(:)   ! LIST OF VARIABLE INFOs
  INTEGER                               :: len
  CHARACTER(GRD_MAXSTRLEN)              :: infile    ! filename

  ! INIT
  IF (PRESENT(file)) THEN
     infile = TRIM(file)
  ELSE
     infile = ''
  END IF
  NULLIFY(kvar)

  IF (ASSOCIATED(ivar)) DEALLOCATE(ivar); NULLIFY(ivar)
  IF (ASSOCIATED(ovar)) DEALLOCATE(ovar); NULLIFY(ovar)
  IF (ASSOCIATED(scl))  DEALLOCATE(scl);  NULLIFY(scl)
  IF (ASSOCIATED(RGT))  DEALLOCATE(RGT);  NULLIFY(RGT)
  IF (ASSOCIATED(RGTstr)) DEALLOCATE(RGTstr); NULLIFY(RGTstr)

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
        ALLOCATE(RGT(nvar), STAT=status)
        CALL ERRMSG(substr,status,4)
        ALLOCATE(RGTstr(nvar), STAT=status)
        CALL ERRMSG(substr,status,5)
        !
        ! INIT
        scl(:) = 1.0_dp
        !
        DO i=1, nvar
           ivar(i) = TRIM(kvar(i)%name)
           ovar(i) = TRIM(kvar(i)%name)
           ! TRY TO GET RGTstr, RGT from attributes
           RGTstr(i) = '...'
           DO j=1, kvar(i)%natts
              IF (TRIM(kvar(i)%att(j)%name) == 'RG_TYPE') THEN
                 RGTstr(i) = string(kvar(i)%att(j)%dat%vc)
              END IF
           END DO
        END DO
        ! CLEAN UP
        DO i=1, nvar
           CALL INIT_NCVAR(kvar(i))
        END DO
        DEALLOCATE(kvar, STAT=status)
        CALL ERRMSG(substr,status,6)
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
     CALL ERRMSG(substr,status,7)
     ALLOCATE(ovar(nvar), STAT=status)
     CALL ERRMSG(substr,status,8)
     ALLOCATE(scl(nvar), STAT=status)
     CALL ERRMSG(substr,status,9)
     ALLOCATE(RGT(nvar), STAT=status)
     CALL ERRMSG(substr,status,10)
     ALLOCATE(RGTstr(nvar), STAT=status)
     CALL ERRMSG(substr,status,11)
     ALLOCATE(lvar(nvar), STAT=status)
     CALL ERRMSG(substr,status,12)
     !
     ! INIT
     scl(:) = 1.0_dp
     !
     ! PARSING
     idx1 = 1
     idx2 = INDEX(TRIM(var), ';')
     DO i=1, nvar-1
        lvar(i) = var(idx1:idx2-1)
        idx1 = idx2+1
        idx2 = idx2+INDEX(var(idx1:),';')
     END DO
     IF (idx2 == (idx1-1)) idx2 = len+1    ! last ';' optional !
     lvar(nvar) = var(idx1:idx2-1)
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
        ! CHECK for destination variable name
        idxl = 1
        idxr = MIN(idx1, idx2, idx3)-1
        ovar(i) = lvar(i)(idxl:idxr)
        ovar(i) = ADJUSTL(TRIM(ovar(i)))
        ! CHECK for source variable name
        idxl = idx1+1
        idxr = MIN(idx2, idx3)-1
        IF (idxl > len) idxl = 1
        IF (idxr > len) idxr = len
        ivar(i) = lvar(i)(idxl:idxr)
        ivar(i) = ADJUSTL(TRIM(ivar(i)))
        ! CHECK for scaling
        idxl = idx2+1
        idxr = idx3-1
        IF (idxr < idxl) idxr = len
        IF (idxl > len) THEN
           scl(i) = 1.0_dp
        ELSE
           READ(lvar(i)(idxl:idxr),*) scl(i)
        END IF
        ! CHECK for RE-GRIDDING TYPE
        idxl = idx3+1
        idxr = idx2-1
        IF (idxr < idxl) idxr = len
        IF (idxl > len) THEN        ! RG_TYPE NOT specified
           RGTstr(i) = '...'
           IF (TRIM(infile) /= '') THEN  ! filename OK
              ALLOCATE(kvar(1), STAT=status)
              CALL ERRMSG(substr,status,13)
              CALL IMPORT_NCVAR(kvar(1), 1, &
                   varname=TRIM(ivar(i)), file=TRIM(infile))
              CALL INIT_NARRAY(kvar(1)%dat)
              DO j=1, kvar(1)%natts
                 IF (TRIM(kvar(1)%att(j)%name) == 'RG_TYPE') THEN
                    RGTstr(i) = string(kvar(1)%att(j)%dat%vc)
                 END IF
              END DO
              ! CLEAN UP
              CALL INIT_NCVAR(kvar(1))
              DEALLOCATE(kvar, STAT=status)
              CALL ERRMSG(substr,status,14)
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
           RGTstr(i) = lvar(i)(idxl:idxr)
           RGTstr(i) = ADJUSTL(TRIM(RGTstr(i)))
        END IF                    ! RG_TYPE (NOT) specified
     END DO  ! LOOP OVER ALL VARIABLES

     ! CLEAN
     IF (ALLOCATED(lvar)) THEN
        DEALLOCATE(lvar, STAT=status)
        CALL ERRMSG(substr,status,15)
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

#ifdef PNCREGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE DISTRIBUTE_VARS_ON_PES(ivar, ovar, scl, RGT, RGTstr)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: ivar(:)   ! input variable names
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: ovar(:)   ! output variable names
  REAL(dp)                    , POINTER :: scl(:)    ! scaling
  INTEGER                     , POINTER :: RGT(:)    ! regridding type
  CHARACTER (LEN=3)           , POINTER :: RGTstr(:) ! regridding type string

  ! LOCAL
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: zivar(:)   ! input variable names
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: zovar(:)   ! output variable names
  REAL(dp)                    , POINTER :: zscl(:)    ! scaling
  INTEGER                     , POINTER :: zRGT(:)    ! regridding type
  CHARACTER (LEN=3)           , POINTER :: zRGTstr(:) ! regridding type string
  !
  INTEGER :: nv, i
  INTEGER :: pe_inc, active_proc
  INTEGER :: n, no
  INTEGER, DIMENSION(:), ALLOCATABLE :: znv

  ! SAVE INPUT
  nv = SIZE(ivar)
  ALLOCATE(zivar(nv))
  ALLOCATE(zovar(nv))
  ALLOCATE(zscl(nv))
  ALLOCATE(zRGT(nv))
  ALLOCATE(zRGTstr(nv))
  !
  zivar(:) = ivar(:)
  zovar(:) = ovar(:)
  zscl(:)  = scl(:)
  zRGT(:)  = RGT(:)
  zRGTstr(:) = RGTstr(:)
  
  ! RE-SET OUTPUT
  DEALLOCATE(ivar);   NULLIFY(ivar)
  DEALLOCATE(ovar);   NULLIFY(ovar)
  DEALLOCATE(scl);    NULLIFY(scl)
  DEALLOCATE(RGT);    NULLIFY(RGT)
  DEALLOCATE(RGTstr); NULLIFY(RGTstr)

  ! determine working PEs
  i_am_worker = .FALSE.
  nproc_work = MIN(nproc, nv, MAX_working_PEs)
  nproc_work = MAX(nproc_work,1)
  IF (ALLOCATED(pe_list)) DEALLOCATE(pe_list)
  ALLOCATE(pe_list(0:nproc_work-1))
  pe_inc = nproc/nproc_work
  DO i=0,nproc_work-1
     pe_list(i) = 0+i*pe_inc
  END DO

  ! distribute variables across working PEs
  ALLOCATE(znv(0:nproc_work-1))
  znv(:) = 0
  active_proc = 0
  DO i=1, nv
     znv(active_proc) = znv(active_proc) + 1
     i_am_worker = i_am_worker .OR. (pe_list(active_proc) == me)
     active_proc = active_proc+1
     active_proc = MOD(active_proc,nproc_work)
  END DO

  no = 0
  DO i=0, nproc_work-1
     n = znv(i)
     IF (pe_list(i) == me) THEN
        ALLOCATE(ivar(n))
        ivar(:) = zivar(no+1:no+n)
        ALLOCATE(ovar(n))
        ovar(:) = zovar(no+1:no+n)
        ALLOCATE(scl(n))
        scl(:) = zscl(no+1:no+n)
        ALLOCATE(RGT(n))
        RGT(:) = zRGT(no+1:no+n)
        ALLOCATE(RGTstr(n))
        RGTstr(:) = zRGTstr(no+1:no+n)
     END IF
     no = no + n
  END DO

  ! CLEAN MEMORY
  DEALLOCATE(zivar);   NULLIFY(zivar)
  DEALLOCATE(zovar);   NULLIFY(zovar)
  DEALLOCATE(zscl);    NULLIFY(zscl)
  DEALLOCATE(zRGT);    NULLIFY(zRGT)
  DEALLOCATE(zRGTstr); NULLIFY(zRGTstr)

END SUBROUTINE DISTRIBUTE_VARS_ON_PES
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE EXPAND_FILENAME(outfile, nvars, me)

  IMPLICIT NONE

  ! I/O
  CHARACTER(len=*),INTENT(INOUT)  :: outfile  
  INTEGER,         INTENT(IN)     :: nvars
  INTEGER,         INTENT(IN)     :: me

  ! LOCAL
  CHARACTER(len=4)         :: cpe

  IF(nproc_work > 1 .AND. outfile /= "") THEN
     WRITE(cpe,'(i3.3,''_'')') me
     outfile = cpe//TRIM(outfile)
  END IF

END SUBROUTINE EXPAND_FILENAME
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE INIT_PNCREGRID

  IMPLICIT NONE

  work_loop = 0
  IF (ALLOCATED(pe_list)) DEALLOCATE(pe_list)
  ALLOCATE(pe_list(0:0))
  pe_list(0) = 0

END SUBROUTINE INIT_PNCREGRID
! ------------------------------------------------------------------

#endif

! ******************************************************************
END MODULE MESSY_NCREGRID_CONTROL
! ******************************************************************
