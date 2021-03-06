#include "messy_main_ppd_bi.inc"

!*************************************************************************
MODULE messy_main_import_grid_bi
!*************************************************************************

  !*************************************************************************
  !  IMPORT GRID
  !   FOR GLOBAL MODELS NCREGRID functionality adopted from import_rgt
  !   SCRIP implemented for data input on global and regional grids
  !
  !  Author: Astrid Kerkweg, Uni Mainz, 2012/2013
  !*************************************************************************
  !  Expansion of IMPORT_GRID to patches (ICON model)
  !  Author: Astrid Kerwkeg, Uni Bonn, 2017
  !*************************************************************************

#if !defined(MBM_QBO)

  USE messy_main_import_grid_tools_bi, ONLY: T_RGTEVENT, T_RGTINIT       &
                                           , RGTMAXACTSTR                &
                                           , RGREAD_NCVAR, RGREAD_NCFILE
  USE messy_main_blather_bi,           ONLY: START_MESSAGE_BI, END_MESSAGE_BI &
                                           , ERROR_BI
  USE messy_main_tools,                ONLY: PTR_4D_ARRAY
  USE messy_main_constants_mem,        ONLY: STRLEN_MEDIUM
  USE messy_main_import_grid,          ONLY: dp, submodstr
  USE messy_main_grid_trafo,           ONLY: GTRFCHAR &
                                           , GTRF_NONE, GTRF_NRGD
  USE messy_main_grid,                 ONLY: t_geohybgrid
  USE messy_main_import,               ONLY: modstr
  USE messy_main_channel,              ONLY: t_chaobj_cpl, STRLEN_OBJECT
  USE messy_main_grid_netcdf,          ONLY: t_multinc
#if ! ( defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE)))
  USE messy_main_grid_trafo,           ONLY: GTRF_SCRP
#endif

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTRINSIC :: ASSOCIATED, LEN_TRIM, NULL, REAL, SIZE, TRIM

  TYPE t_import_set
     CHARACTER(LEN=RGTMAXACTSTR)     :: nml = '' ! file with ONE regrid-namelist
     CHARACTER(LEN=RGTMAXACTSTR)     :: var = '' ! variable(tracer) to select
     CHARACTER(LEN=RGTMAXACTSTR)     :: file= '' ! file to select
     ! string containing heights for Nx2D
     CHARACTER(LEN=RGTMAXACTSTR)     :: heights=''
     REAL(DP), DIMENSION(:), POINTER :: z => NULL()! emission heights ...
                                               ! ... for multi level emissions
     CHARACTER(LEN=RGTMAXACTSTR)     :: presslevs=''
     REAL(DP), DIMENSION(:), POINTER :: p => NULL()! emission press. levels ...
                                               ! ... for multi level emissions
     INTEGER                         :: RGREAD = 0 ! RGREAD_NCVAR, RGREAD_NCFILE
     INTEGER                         :: nvar   = 0 ! number of variables
     ! EMISSION SET OK?
     LOGICAL                         :: lok   = .false.
     ! EMISSION TYPE (1: 2D surface, 2: 3D volume, 3: Nx2D multi level...
     !                4: index regridding, 5: 4D )
     INTEGER                         :: itype = 0
     ! NUMBER OF VALID EMISSION FIELDS (=CHANNEL OBJECTS) IN SET
     INTEGER                         :: nf    = 0
     ! RANGE OF INDEX LEVELS FOR INDEX REGRIDDING "min:max"
     INTEGER,                POINTER :: IXF(:,:) => NULL()
     ! POINTER TO CHANNEL OBJECT FOR INPUT DATA
     TYPE(PTR_4D_ARRAY),  DIMENSION(:), POINTER :: field => NULL()
     ! name of channel objects (required for restart check)
     TYPE(t_chaobj_cpl),  DIMENSION(:), POINTER :: fieldname => NULL()
     ! INTERPOL METHOD: NONE, NCREGRID, SCRIP
#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
     ! set default to NCREGRID to keep functionality of old import_rgt namelists
     INTEGER                     :: GTRF_TYPE = GTRF_NRGD
#else
     ! set default to SCRIP to be functional
     INTEGER                     :: GTRF_TYPE = GTRF_SCRP
#endif
     ! TARGET GRID NAME
     CHARACTER(LEN=STRLEN_MEDIUM):: ogrid_name = ''
     ! GRID IDs
     INTEGER                     :: ogrid_ID = -99 ! output (target) grid
     INTEGER                     :: igrid_ID = -99 ! input grid
     ! DATA ID (REQUIRED IN CASE OF SCRIP interpolation)
     INTEGER                     :: SCRIP_DATA_ID = -99
     !
     ! STORE INFO OF multi-netcdf descriptpr
     TYPE(t_multinc)             :: mnc
     !
  END TYPE t_import_set

  TYPE(t_rgtinit),    DIMENSION(:), POINTER :: RGTINIT     => NULL()
  INTEGER                                   :: NIMPINIT_SET

  TYPE(t_rgtevent),   DIMENSION(:), POINTER :: RGT         => NULL()
  INTEGER                                   :: NIMPORT_SET

  TYPE(t_import_set), DIMENSION(:), POINTER :: IMPORT_SET  => NULL()
  INTEGER                                   :: NIMPALL_SET
  !
  ! Skip code parts of not required
  ! is 3D data (vertical interpolation) required?
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: L_req3Dipol

  PUBLIC :: import_grid_initialize
  PUBLIC :: import_grid_init_memory
  PUBLIC :: import_grid_init_tracer
  PUBLIC :: import_grid_global_start
  PUBLIC :: import_grid_free_memory
  !PRIVATE :: add_var_attributes

CONTAINS

  ! ------------------------------------------------------------------------
  SUBROUTINE  IMPORT_GRID_INITIALIZE

    ! import_grid MODULE ROUTINE (BASE MODEL INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    ! INITIALIZATION OF import_grid SPECIFIC EVENTS FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, Mar 2004

    ! ECHAM5/MESSy
    USE MESSY_MAIN_MPI_BI,               ONLY: p_parallel_io, p_pe
    USE MESSY_MAIN_IMPORT_GRID_TOOLS_BI, ONLY: RGTEVENT_INIT_NML    &
                                             , RGTINIT_INIT_NML     &
                                             , RGTMAXACTSTR
    USE MESSY_MAIN_TIMER_BI,             ONLY: main_timer_set_domain
    USE MESSY_MAIN_CHANNEL_BI,           ONLY: n_dom

    USE MESSY_MAIN_IMPORT_GRID,          ONLY: PARSE_STR
    USE MESSY_MAIN_GRID_NETCDF,          ONLY: my_rank
    USE MESSY_MAIN_TOOLS,                ONLY: MISSING

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'import_grid_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: i
    INTEGER                     :: idx
    CHARACTER(LEN=RGTMAXACTSTR) :: act = ''
    CALL START_MESSAGE_BI(submodstr, 'GLOBAL SETUP',substr)

    status = 0

    ALLOCATE(L_req3Dipol(n_dom))
    L_req3Dipol(:)= .FALSE.

    ! INITIALIZE DEFAULT RGT-INIT
    CALL RGTINIT_INIT_NML(RGTINIT, modstr)
    IF (status /= 0 .AND. status /= MISSING) &
         CALL error_bi ('ERROR IN RGTINIT_INIT_NML',substr)

    NIMPINIT_SET = SIZE(RGTINIT)

    ! INITIALIZE DEFAULT RGT-EVENTS
    CALL RGTEVENT_INIT_NML(RGT, modstr)
    IF (status /=0) CALL error_bi ('ERROR IN RGTEVENT_INIT_NML',substr)

    ! MEMORY SPACE FOR IMPORT-HANDLING
    NIMPORT_SET = SIZE(RGT)

    NIMPALL_SET = NIMPORT_SET+NIMPINIT_SET
    ALLOCATE(IMPORT_SET(NIMPALL_SET))

    ! PARSE ACTION STRING
    IF (p_parallel_io) THEN
       WRITE(*,*) '---------------------------------------------------------'
       WRITE(*,*) 'NUMBER OF IMPORT-SETS event :', NIMPORT_SET
       WRITE(*,*) 'NUMBER OF IMPORT-SETS init  :', NIMPINIT_SET
       WRITE(*,*) 'NUMBER OF IMPORT-SETS all   :', NIMPALL_SET
       WRITE(*,*) '-------------------------------------------------------'
    END IF
    DO i=1, NIMPALL_SET
       IF (i<=NIMPORT_SET) THEN
          act = RGT(i)%io%act
       ELSE
          act = RGTINIT(i-NIMPORT_SET)%act
       END IF
       CALL PARSE_STR(status, RGTMAXACTSTR, act        &
            , IMPORT_SET(i)%nml ,IMPORT_SET(i)%var     &
            , IMPORT_SET(i)%file                       &
            , IMPORT_SET(i)%z, IMPORT_SET(i)%heights   &
            , IMPORT_SET(i)%p, IMPORT_SET(i)%presslevs &
            , IMPORT_SET(i)%GTRF_TYPE                  &
            , IMPORT_SET(i)%ogrid_name)
       ! DEFAULT NAMELIST
       IF (TRIM(IMPORT_SET(i)%nml) ==  '') &
            IMPORT_SET(i)%nml = modstr//'.nml'
       idx = LEN_TRIM(IMPORT_SET(i)%nml)
       ! remove '.nml' suffix since it is appended in NCREGRID tools
       IF (IMPORT_SET(i)%nml(idx-3:idx) == ".nml") &
            IMPORT_SET(i)%nml = IMPORT_SET(i)%nml(1:idx-4)

       ! CHECK CONSISTENCY
       ! (1) FILE AND VAR -> NCVAR
       !              VAR -> NCVAR
       !     FILE         -> NCFILE
       !                  -> NCFILE
       !
       IF (TRIM(IMPORT_SET(i)%var) /= '') THEN
          ! PROCESS ONE VARIABLE
          IMPORT_SET(i)%RGREAD = RGREAD_NCVAR
          IF (TRIM(IMPORT_SET(i)%file) /= '') THEN
             IF (p_parallel_io) THEN
                WRITE(*,*) substr//' - WARNING:'
                WRITE(*,*) '   VAR and FILE set IN RGT-EVENT'//&
                     &' ACTION STRING ',i
                WRITE(*,*) '   ( '//TRIM(RGT(i)%io%act)//' )'
                WRITE(*,*) '   -> FILE SPECIFIER IGNORED !!!'
             END IF
             IMPORT_SET(i)%file = ''
          END IF
       ELSE
          ! PROCESS ONE NC-FILE
          IMPORT_SET(i)%RGREAD = RGREAD_NCFILE
       END IF

       ! OUTPUT RESULT
       IF (p_parallel_io) THEN
          WRITE(*,*) 'IMPORT-SET NUMBER  : ', i
#ifdef ICON
          IF (i<= NIMPORT_SET) THEN
             WRITE(*,*) 'IMPORT-SET DOMAIN  : ', RGT(i)%domain_idx
          ELSE
             WRITE(*,*) 'IMPORT-SET DOMAIN  : ', &
                  RGTINIT(i-NIMPORT_SET)%domain_idx
          END IF
#endif
          WRITE(*,*) 'IMPORT-SET ACTION  : ', act
          WRITE(*,*) 'IMPORT-SET NML     : ' &
               , TRIM(IMPORT_SET(i)%nml)//'.nml'
          WRITE(*,*) 'IMPORT-SET VAR     : ', TRIM(IMPORT_SET(i)%var)
          WRITE(*,*) 'IMPORT-SET FILE    : ', TRIM(IMPORT_SET(i)%file)
          SELECT CASE(IMPORT_SET(i)%RGREAD)
          CASE(RGREAD_NCFILE)
             WRITE(*,*) 'IMPORT_SET RGREAD  : NCFILE'
          CASE(RGREAD_NCVAR)
             WRITE(*,*) 'IMPORT_SET RGREAD  : NCVAR'
          CASE DEFAULT
             CALL ERROR_BI('UNKNOWN RGREAD METHOD',substr)
          END SELECT
          IF (ASSOCIATED(IMPORT_SET(i)%z)) THEN
             WRITE(*,*) 'IMPORT LEVELS [m] : ',IMPORT_SET(i)%z
          END IF
          IF (ASSOCIATED(IMPORT_SET(i)%p)) THEN
             WRITE(*,*) 'IMPORT PRESSURE LEVELS [Pa] : ' &
                  ,IMPORT_SET(i)%p
          END IF
#if !defined(ECHAM5) && !defined(CESM1)
          IF (IMPORT_SET(i)%GTRF_TYPE == GTRF_NRGD) &
               CALL ERROR_BI('NCREGRID interpolation only possible for ECHAM5' &
               , substr)
#endif
          WRITE(*,*) 'IMPORT-SET INTERPOLATION  : ' &
               , GTRFCHAR(IMPORT_SET(i)%GTRF_TYPE)

          IF (LEN_TRIM(IMPORT_SET(i)%ogrid_name) > 0) &
               WRITE(*,*) 'IMPORT-SET TARGET GRID   : ' &
               ,TRIM(IMPORT_SET(i)%ogrid_name)

          WRITE(*,*) '--------------------------------------------'
       END IF

       SELECT CASE (status)
       CASE(0)
          ! PROCESS NAMELISTS in init_memory,
       CASE(1)
          CALL ERROR_BI('ERROR IN RGT-EVENT ACTION-STRING',substr)
       CASE(2)
          CALL ERROR_BI('EMPTY SPECIFICATION IN RGT-EVENT '//&
               &'ACTION STRING',substr)
       CASE(3)
          CALL ERROR_BI('UNKNOWN SPECIFIER IN RGT-EVENT '//&
               &'ACTION-STRING',substr)
       CASE(6)
          CALL ERROR_BI('ERROR IN Z-SPECIFIER IN RGT-EVENT '//&
               &'ACTION-STRING',substr)
       CASE DEFAULT
          CALL ERROR_BI('ERROR IN RGT-EVENT ACTION-STRING',substr)
       END SELECT
    END DO

    my_rank = p_pe
    CALL main_timer_set_domain(1)

    CALL END_MESSAGE_BI(submodstr, 'GLOBAL SETUP',substr)

  END SUBROUTINE IMPORT_GRID_INITIALIZE
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE IMPORT_GRID_INIT_MEMORY

#if defined(ECHAM5)
    USE messy_main_mpi_bi,        ONLY: P_BCAST, p_io
#endif
    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    USE messy_main_grid_def_mem_bi,      ONLY: nlev, nproma &
#ifndef ICON
                                             , ngpblks       &
#endif
                                             , BASEGRID_ID
    USE messy_main_import_grid_tools_bi, ONLY: RGTOOL_BI_READ_NCVAR &
                                       , RGTOOL_BI_READ_NCFILE      &
                                       , RGTEVENT_TO_CHANNEL
    USE messy_main_channel_error_bi,     ONLY: CHANNEL_HALT
#ifndef ICON
    USE messy_main_channel_bi,    ONLY: GP_3D_MID, GP_3D_1LEV   &
                                      , gp_cnt, gp_meml, gp_memu
#endif
    USE messy_main_channel_bi,    ONLY: gp_start  &
                                      , gp_nseg &
                                      , DIMID_LON, DIMID_LAT, DC_GP, DC_BC &
                                      , DIMID_LEV
#ifdef ICON
    USE messy_main_channel_bi,    ONLY: DIMID_ONE, DIMID_NCELLS        &
                                      , UCELL_3D_1LEV, UCELL_3D_MID
    USE messy_main_bmluse_bi,     ONLY: p_patch
#endif
    USE messy_main_grid_netcdf,   ONLY: GRD_MAXSTRLEN
    USE messy_main_grid,          ONLY: GRID_ERROR, LOCATE_GEOHYBGRID
    USE messy_main_channel,       ONLY: NEW_CHANNEL, NEW_CHANNEL_OBJECT &
                                      , NEW_ATTRIBUTE
    USE messy_main_tools,         ONLY: INT2STR, STRCRACK, STR2NUM
    USE messy_main_channel_dimensions, ONLY: NEW_DIMENSION              &
                                           , GET_DIMENSION_INFO         &
                                           , ADD_DIMENSION_VARIABLE     &
                                           , ADD_DIMENSION_VARIABLE_ATT &
                                           , DIMID_UNDEF                &
                                           , GET_DIMVAR_INFO

    USE messy_main_channel_repr,       ONLY: GET_REPRESENTATION_INFO   &
                                           , GET_REPRESENTATION        &
                                           , NEW_REPRESENTATION        &
                                           , AUTO, IRANK
#ifndef ICON
    USE messy_main_channel_repr,       ONLY: SET_REPRESENTATION_DECOMP &
                                           , PIOTYPE_COL  &
                                           , REPR_DEF_AXES
#endif
    USE messy_main_channel_attributes, ONLY: t_att_ptr
    USE messy_main_channel_mem,        ONLY: dom_curid
    USE messy_main_constants_mem,      ONLY: STRLEN_ULONG, STRLEN_MEDIUM

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER      :: substr = 'import_grid_init_memory'
    INTEGER                          :: status
    INTEGER                          :: i, jv
    REAL(dp), DIMENSION(:,:,:,:,:),             POINTER :: dat => NULL()
    CHARACTER(LEN=GRD_MAXSTRLEN), DIMENSION(:), POINTER :: var => NULL()
    LOGICAL                          :: lok
    INTEGER                          :: mvar     ! number ov variables
    INTEGER                          :: mlev     ! number of levels in IMPORT
                                                 ! channel object (2D or 3D)
    INTEGER                          :: ip       ! length of parameter dim.
    INTEGER                          :: reprid   ! REPRESENTATION ID
    INTEGER                          :: DIMID_VERT
    INTEGER                          :: DIMID_PAR
    INTEGER                          :: nseg    = 0  ! max no. of segments
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()
    CHARACTER(LEN=*), PARAMETER      :: ch_str = 'import_grid'
    CHARACTER(LEN=3)                 :: numstr
    INTEGER                          :: xdim, ydim, zdim, ndim
    INTEGER                          :: xdim_ID, ydim_ID, zdim_ID, ndim_ID
    CHARACTER(LEN=4)                 :: dimstr
    CHARACTER(LEN=21)                :: repr_str
    LOGICAL                          :: lnumdim = .FALSE.
    TYPE(t_att_ptr), DIMENSION(:), POINTER :: att => NULL()
    INTEGER                          :: xnpar, xnlev, xnlon, xnlat
    INTEGER                          :: np
    CHARACTER(LEN=2)                 :: spid
    LOGICAL                          :: lrestart = .FALSE.
    INTEGER                                            :: curcnt = 0
    INTEGER                                            :: curixf = 0
    INTEGER                                            :: reprcnt = 1000
    INTEGER                                            :: ic, l
    CHARACTER(LEN=3)                                   :: cntstr  = ' '
    CHARACTER(LEN=STRLEN_MEDIUM)                       :: dimname = ' '
    CHARACTER(LEN=STRLEN_ULONG)                        :: levunit = ' '
    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER :: chlevs  => NULL()
    REAL(dp),                    DIMENSION(:), POINTER :: levs    => NULL()
    REAL(dp),                    DIMENSION(:), POINTER :: dimlevs => NULL()
    INTEGER                                            :: idxrun  = 0
    LOGICAL                                            :: l_mkdim = .FALSE.
    LOGICAL                                            :: l_cyc   = .FALSE.
    INTEGER                                            :: dlen    = 0
    CHARACTER(LEN=3)                                   :: idxstr  = ' '
#if ! defined(ICON)
    CHARACTER(LEN=*), PARAMETER                        :: reprstr = 'GP_'
#else
    CHARACTER(LEN=*), PARAMETER                        :: reprstr = 'UCELL_'
#endif

    INTEGER                                            :: dom_idx
    CHARACTER(LEN=STRLEN_OBJECT)                       :: name
    LOGICAL                                            :: lstatic


    CALL START_MESSAGE_BI(submodstr, 'MEMORY SETUP',substr)

    np  = dom_curid

    CALL int2str(spid,np,'0')
    curcnt  = 0
    reprcnt = 0

    ! define new channel
    CALL NEW_CHANNEL(status, submodstr)
    CALL CHANNEL_HALT(substr, status)

    ! SAVE EVENT / STEPPERS TO CHANNEL FOR RESTART
    CALL RGTEVENT_TO_CHANNEL(submodstr, rgt)

    import_set_loop: DO i=1, NIMPALL_SET

       IF (i<=NIMPORT_SET) THEN
          dom_idx = RGT(i)%domain_idx
          name    = RGT(i)%io%cnt%name
          lstatic = .FALSE.
       ELSE
          dom_idx = RGTINIT(i-NIMPORT_SET)%domain_idx
          name    = RGTINIT(i-NIMPORT_SET)%name
          lstatic = .TRUE.
       END IF

       IF (dom_idx /= dom_curid) CYCLE

       IF (p_parallel_io) THEN
          WRITE(*,*) '========================================================'
          WRITE(*,*) 'IMPORT-SET       : ',i
       END IF

       ! define output grid: DEFAULT: base model grid
       ! TODO: add variable output grid
       IF (IMPORT_SET(i)%GTRF_TYPE == GTRF_NONE) THEN
          IMPORT_SET(i)%ogrid_ID = -99
          !NULLIFY(outarea)
       ELSE IF (LEN_TRIM(IMPORT_SET(i)%ogrid_name) > 0) THEN
          CALL LOCATE_GEOHYBGRID(status &
               , TRIM(IMPORT_SET(i)%ogrid_name), ID=IMPORT_SET(i)%ogrid_ID)
          IF (status /= 0) CALL ERROR_BI(grid_error(status), substr)
       ELSE
          IMPORT_SET(i)%ogrid_ID = BASEGRID_ID(np)
       ENDIF

       SELECT CASE(IMPORT_SET(i)%RGREAD)
       CASE(RGREAD_NCFILE)
          !
          CALL RGTOOL_BI_READ_NCFILE( TRIM(IMPORT_SET(i)%nml)   &
               , TRIM(IMPORT_SET(i)%file), 1, dat               &
               , IMPORT_SET(i)%GTRF_TYPE                        &
               , IMPORT_SET(i)%IXF                              &
               , IMPORT_SET(i)%ogrid_ID, IMPORT_SET(i)%igrid_ID &
               , var, SDID=IMPORT_SET(i)%SCRIP_DATA_ID          &
               , lrg = .false. , lok = lok                      &
               , att=att                                        &
               , ldatreq = .FALSE.                              &
               , nlon=xnlon, nlat=xnlat                         &
               , nlev=xnlev, npar=xnpar                         &
               , mnc=IMPORT_SET(i)%mnc )
          !
       CASE(RGREAD_NCVAR)
          !
          CALL RGTOOL_BI_READ_NCVAR( TRIM(IMPORT_SET(i)%nml)    &
               , TRIM(IMPORT_SET(i)%var), 1, dat                &
               , IMPORT_SET(i)%GTRF_TYPE                        &
               , IMPORT_SET(i)%IXF                              &
               , IMPORT_SET(i)%ogrid_ID, IMPORT_SET(i)%igrid_ID &
               , SDID=IMPORT_SET(i)%SCRIP_DATA_ID               &
               , lrg = .false. , lok = lok                      &
               , att=att                                        &
               , ldatreq = .FALSE.                              &
               , nlon=xnlon, nlat=xnlat                         &
               , nlev=xnlev, npar=xnpar                         &
               , mnc=IMPORT_SET(i)%mnc )
          ALLOCATE(var(1))
          var(1) = TRIM(IMPORT_SET(i)%var)
          !
       CASE DEFAULT
          !
          ! SHOULD BE NEVER REACHED
          CALL ERROR_BI('UNKNOWN RGREAD METHOD',substr)
          !
       END SELECT

#if defined(ECHAM5)
       CALL P_BCAST(IMPORT_SET(i)%igrid_id, p_io)
       CALL P_BCAST(IMPORT_SET(i)%ogrid_id, p_io)
#endif

       IMPORT_SET(i)%lok = lok

       IF (.NOT. lok) CALL ERROR_BI(' NOT OK ', substr)

       IF (lok) THEN
          ! NUMBER OF VARIABLES
          mvar = SIZE(var)
          IMPORT_SET(i)%nvar = mvar
          IF (p_parallel_io) THEN
             WRITE(*,*) 'IMPORT-SET       : ',i, 'patch:'//spid
             WRITE(*,*) '  NO. OF VARIABLES: ',mvar, 'patch:'//spid
          END IF
          ! ALLOCATE POINTER TO CHANNEL OBJECTS FOR INPUT
          ALLOCATE(IMPORT_SET(i)%field(mvar))
          ALLOCATE(IMPORT_SET(i)%fieldname(mvar))

          DO jv=1, mvar
             NULLIFY(IMPORT_SET(i)%field(jv)%ptr)
             IMPORT_SET(i)%fieldname(jv)%obj = ' '
             IMPORT_SET(i)%fieldname(jv)%cha = ch_str
          END DO

          ! ************************
          ! DETECT TYPE OF IMPORT:
          ! ************************
          !
          ! NCVAR : dat(#lon, #lev, #param, #lat, 1)
          ! NCFILE: dat(#lon, #lev, #param, #lat, #nvar)
          !
          !        SIZE    SIZE     SIZE     ASSOCIATED
          !        (var)   (dat,2)  (dat,3)  (z)
          !                = mlev   = ...
          !        ======================================================
          ! NCVAR    1       1       1       no          -> 2D
          ! NCVAR    1      nlev     1       no          -> 3D
          ! NFILE   1,nvar   1      1,nvar   no          -> 2D
          ! NFILE   1,nvar  nlev    1,nvar   no          -> 3D
          ! NCVAR    1       1      1,N      no          ->
          ! NCVAR    1       1      1,N      yes         -> Nx2D
          !
          !
          mlev = xnlev
          ip   = xnpar
          lrestart = .FALSE.

          level_analysis: IF (IMPORT_SET(i)%GTRF_TYPE == GTRF_NONE) THEN
             IMPORT_SET(i)%itype = 0 ! import raw data
          ELSE IF (mlev > 1) THEN
             !
             IF (ASSOCIATED(IMPORT_SET(i)%p)) THEN
                IMPORT_SET(i)%itype = 3 ! multilayer emissions pressure levels
             ELSE IF (ip == 1) THEN
                ! IF THE NUMBER OF LEVELS IS > 1, a FIELD in
                ! GP_3D_MID  REPRESENTATION IS ASSUMED,
                !
                IMPORT_SET(i)%itype = 2       ! 3D
                lrestart = .TRUE.
                 L_req3Dipol = .TRUE.
             ELSE
                IMPORT_SET(i)%itype = 5       ! 4D
                lrestart = .TRUE.
             ENDIF
          ELSE ! (mlev == 1)
             ! IF THE NUMBER OF LEVELS IS 1, EITHER A HORIZONTAL-2D-FIELD OR A
             ! MULTI LEVEL 2D-FIELD (Nx2D) IS ASSUMED.
             ! NO VERTICAL AXIS HAS BEEN SPECIFIED IN THE RESPECTIVE
             ! &REGRID NAMELIST. THE LEVELS FOR Nx2D FIELDS ARE
             ! THEN TREATED BY NCREGRID AS 'INVARIANT AXIS' AND
             ! OCCUR AT RANK 3 OF THE DATA FIELD dat.
             ! HOWEVER, SINCE NCFILE PROCESSING CAN ALSO BE APPLIED TO
             ! netCDF FILES WITH ONLY ONE VARIABLE (RESULTING IN
             ! nvar=1 AND SIZE(dat,3)=1), AND FURTHER Nx2D VARIABLES CAN
             ! HAVE AN ARBITRARY NUMBER OF LEVELS (INCLUDING 1 !),
             ! THE ONLY POSSIBLE CRITERIUM TO TRIGGER Nx2D HANDLING
             ! IS THE PRESENCE OF A HEIGHT SPECIFIER (Z) IN THE NAMELIST ...
             !
             IF (ASSOCIATED(IMPORT_SET(i)%z)) THEN
                ! Nx2D (multi level field)
                IMPORT_SET(i)%itype = 3  ! Nx2D (MULTI LEVEL FIELD)
                !
                ! CHECK CONSISTENCY (number of levels in namelist)
                IF (ip /= SIZE(IMPORT_SET(i)%z)) THEN
                   IF (p_parallel_io) THEN
                      WRITE(*,*) substr//' - ERROR:'
                      WRITE(*,*) 'NUMBER OF IMPORT LEVELS IN NAMELIST (',&
                           SIZE(IMPORT_SET(i)%z),')'
                      WRITE(*,*) 'INCONSISTENT WITH DIMENSION LENGTH IN'
                      WRITE(*,*) 'netCDF FILE (',ip,') !!!'
                   END IF
                   CALL ERROR_BI('HEIGHT SPECIFIER INCONSISTENT',substr)
                END IF
             ELSE IF (ASSOCIATED(IMPORT_SET(i)%p)) THEN
                ! Nx2D (multi level field)
                IMPORT_SET(i)%itype = 3  ! Nx2D (MULTI LEVEL FIELD)
                !
                ! CHECK CONSISTENCY (number of levels in namelist)
                IF (ip /= SIZE(IMPORT_SET(i)%p)) THEN
                   IF (p_parallel_io) THEN
                      WRITE(*,*) substr//' - ERROR:'
                      WRITE(*,*) 'NUMBER OF IMPORT LEVELS IN NAMELIST (',&
                           SIZE(IMPORT_SET(i)%p),')'
                      WRITE(*,*) 'INCONSISTENT WITH DIMENSION LENGTH IN'
                      WRITE(*,*) 'netCDF FILE (',ip,') !!!'
                   END IF
                   CALL ERROR_BI('PRESSURE LEVEL SPECIFIER INCONSISTENT',substr)
                END IF
             ELSE IF (ANY(IMPORT_SET(i)%IXF /= 0)) THEN
                ! INDEX REGRIDDING 2D => 3D FIELD
                IMPORT_SET(i)%itype = 4
             ELSE IF (ip > 1) THEN
                IMPORT_SET(i)%itype = 6 ! 4D with mlev = 1, nparam > 1
             ELSE
                ip = 1
                IMPORT_SET(i)%itype = 1        ! 2D FIELD
             END IF
          END IF level_analysis

          SELECT CASE(IMPORT_SET(i)%itype)
          CASE(0)

             ! TODO: definition of lon, lat dimension variables
             !       if accessible from ingrid
             ! IMPORT RAW DATA SET
             !
             ! Determine dimensions
             xdim = xnlon
             ydim = xnlat
             zdim = xnlev
             ndim = xnpar

             IF (ndim > 1 .AND. zdim > 1) &
                  CALL ERROR_BI('raw_field has to many ranks', substr)

             ! 1st dimension
             CALL INT2STR(dimstr,xdim)
             CALL GET_DIMENSION_INFO( status, 'DIM_'//dimstr//'LEN'   &
                                    , id=xdim_id, len=xdim, dom_id=np)
             IF (status /= 0) THEN
                CALL NEW_DIMENSION( status, xdim_id &
                                  , 'DIM_'//dimstr//'LEN', xdim, dom_id=np)
                CALL CHANNEL_HALT(substr,status)
                IF (p_parallel_io) &
                     write (*,*) 'DEFINE DIMENSION ','DIM_'//dimstr//'LEN'
             ELSE
                IF (p_parallel_io) write (*,*) &
                     'DIMENSION ','DIM_'//dimstr//'LEN exists already'
             ENDIF
             repr_str='RAW_'//dimstr

             ! 2nd dimension
             CALL INT2STR(dimstr,ydim)
             CALL GET_DIMENSION_INFO( status, 'DIM_'//dimstr//'LEN'   &
                                    , id=ydim_id, len=ydim, dom_id=np)
             IF (status /= 0) THEN
                CALL NEW_DIMENSION( status, ydim_id &
                                  , 'DIM_'//dimstr//'LEN', ydim, dom_id=np)
                CALL CHANNEL_HALT(substr,status)
                IF (p_parallel_io) &
                     write (*,*) 'DEFINE DIMENSION ','DIM_'//dimstr//'LEN'
             ELSE
                IF (p_parallel_io) write (*,*) &
                     'DIMENSION ','DIM_'//dimstr//'LEN exists already'
             ENDIF
             repr_str=TRIM(repr_str)//dimstr

             IF (ndim == 1) THEN
                lnumdim = .FALSE.
                ! 3rd dimension
                CALL INT2STR(dimstr,zdim)
                CALL GET_DIMENSION_INFO( status, 'DIM_'//dimstr//'LEN'   &
                     , id=zdim_id, len=zdim, dom_id=np)
                IF (status /= 0) THEN
                   CALL NEW_DIMENSION( status, zdim_id &
                        , 'DIM_'//dimstr//'LEN', zdim, dom_id=np)
                   CALL CHANNEL_HALT(substr,status)
                   IF (p_parallel_io) &
                        write (*,*) 'DEFINE DIMENSION ','DIM_'//dimstr//'LEN'
                ELSE
                   IF (p_parallel_io) write (*,*) &
                        'DIMENSION ','DIM_'//dimstr//'LEN exists already'
                ENDIF
                repr_str=TRIM(repr_str)//dimstr//'Z'
             ELSE
                ! 4th dimension
                lnumdim = .TRUE.
                CALL INT2STR(dimstr,ndim)
                CALL GET_DIMENSION_INFO( status, 'DIM_'//dimstr//'LEN'   &
                     , id=ndim_id, len=ndim, dom_id=np)
                IF (status /= 0) THEN
                   CALL NEW_DIMENSION( status, ndim_id &
                        , 'DIM_'//dimstr//'LEN', ndim, dom_id=np)
                   CALL CHANNEL_HALT(substr,status)
                   IF (p_parallel_io) &
                        write (*,*) 'DEFINE DIMENSION ','DIM_'//dimstr//'LEN'
                ELSE
                   IF (p_parallel_io) write (*,*) &
                        'DIMENSION ','DIM_'//dimstr//'LEN exists already'
                ENDIF
                repr_str=TRIM(repr_str)//dimstr//'N'
             ENDIF

             ! DEFINE CHANNEL OBJECT REPRESENTATION

             CALL GET_REPRESENTATION_INFO(status, repr_str, id=reprid &
                                         , dom_id=np)

             IF (status == 2003) THEN
                ! in this case the representation was not jet made
                nseg = gp_nseg

                ALLOCATE(start(nseg,IRANK))
                ALLOCATE(cnt(nseg,IRANK))
                ALLOCATE(meml(nseg,IRANK))
                ALLOCATE(memu(nseg,IRANK))

                start(:,:)=1
                cnt(:,:)  =1
                meml(:,:) =1
                memu(:,:) =1

                IF (lnumdim) THEN
                   CALL NEW_REPRESENTATION(status, reprid, repr_str      &
                        , rank = 3, link = 'xxx-', dctype = DC_BC        &
                        , dimension_ids = (/ xdim_id, ydim_id, ndim_id/) &
                        , ldimlen       = (/ xdim, ydim, ndim/)          &
                        , output_order  = (/ 1,2,3 /)                  &
                        , axis          = 'XYN-'                         &
                        , dom_id=np)
                   CALL CHANNEL_HALT(substr//'01', status)
                   start(:,:) = gp_start(:,:)

                   cnt(:,1)  = xdim
                   cnt(:,2)  = ydim
                   cnt(:,4)  = 1
                   cnt(:,3)  = ndim
                   memu(:,1) = xdim
                   memu(:,2) = ydim
                   memu(:,4) = 1
                   memu(:,3) = ndim

                ELSE
                   CALL NEW_REPRESENTATION(status, reprid, repr_str      &
                        , rank = 3, link = 'xxx-', dctype = DC_BC        &
                        , dimension_ids = (/ xdim_id, ydim_id, zdim_id/) &
                        , ldimlen       = (/ xdim, ydim, zdim/)          &
                        , output_order  = (/ 1,2,3 /)                  &
                        , axis          = 'XYZ-'                         &
                        , dom_id=np)
                   CALL CHANNEL_HALT(substr//'02', status)
                   start(:,:) = gp_start(:,:)

                   cnt(:,1)  = xdim
                   cnt(:,2)  = ydim
                   cnt(:,3)  = zdim
                   cnt(:,4)  = 1
                   memu(:,1) = xdim
                   memu(:,2) = ydim
                   memu(:,3) = zdim
                   memu(:,4) = 1
                ENDIF
#ifndef ICON
                CALL SET_REPRESENTATION_DECOMP(status, reprid &
                     , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
                CALL CHANNEL_HALT(substr, status)
#endif
                DEALLOCATE(start)
                DEALLOCATE(cnt)
                DEALLOCATE(meml)
                DEALLOCATE(memu)
             ELSE
                CALL CHANNEL_HALT(substr, status)
             ENDIF

             ! DIAGNOSTIC OUTPUT
             IF (p_parallel_io) THEN
                   WRITE(*,*) '  IMPORT TYPE   : IMPORT RAW FILE'
                   WRITE(*,*) '  DIMENSIONS X / Y / Z  --> ' &
                        , xdim,' / ',ydim,' / ',zdim
             END IF
             ! ----------------------------------------------------------
           CASE(1)
             ! 2D FIELD
             !
             ! NUMBER OF LEVELS IN CHANNEL OBJECT
             mlev = 1
#if defined(ICON)
             reprid = UCELL_3D_1LEV(np)
#else
             reprid = GP_3D_1LEV
#endif

             ! DIAGNOSTIC OUTPUT
             IF (p_parallel_io) THEN
                WRITE(*,*) '  IMPORT TYPE : regular 2D FIELD'
             END IF
             ! ----------------------------------------------------------
          CASE(2)
             ! 3D regular FIELD

             ! NUMBER OF LEVELS IN CHANNEL OBJECT
             mlev = nlev
#if defined(ICON)
             reprid = UCELL_3D_MID(np)
#else
             reprid = GP_3D_MID
#endif

             ! DIAGNOSTIC OUTPUT
             IF (p_parallel_io) THEN
                WRITE(*,*) '  IMPORT TYPE : regular 3D FIELD'
             END IF
             ! ----------------------------------------------------------
          CASE(3)
             ! Nx2D (MULTI LEVEL FIELD)
             !
             ! NUMBER OF LEVELS IN CHANNEL OBJECT
             mlev = ip
             CALL INT2STR(numstr,ip)
             IF (ASSOCIATED(IMPORT_SET(i)%z)) THEN
                dimname='heights'
                CALL STRCRACK(TRIM(ADJUSTL(IMPORT_SET(i)%heights)) &
                     , ',', chlevs, l)
                levunit = 'm'
             ELSE
                dimname='pressurelevels'
                CALL STRCRACK(TRIM(ADJUSTL(IMPORT_SET(i)%presslevs))&
                     , ',', chlevs, l)
                levunit = 'Pa'
             END IF

             l_mkdim = .TRUE.
             idxrun = 0
             ALLOCATE(levs(SIZE(chlevs)))
             DO ic = 1, SIZE(chlevs)
                CALL str2num(chlevs(ic), levs(ic))
             END DO

             DO
                IF (ASSOCIATED(dimlevs)) THEN
                   DEALLOCATE(dimlevs) ; NULLIFY(dimlevs)
                END IF
                idxrun = idxrun +1
                CALL INT2STR(idxstr, idxrun)
                CALL GET_DIMENSION_INFO(status       &
                     , 'DIM_'//numstr//'LEV'//idxstr &
                     , DIMID_VERT, dlen, dom_id=np)
                IF (status == 905) EXIT ! DIMENSION DOES NOT EXIST
                CALL channel_halt(substr,status)

                IF (dlen == ip) THEN
                   !DIMENSION has same length, test attribute
                   CALL get_dimvar_info(status,DIMID_VERT &
                        , dimlevs, firstname=TRIM(dimname))
                   IF (status /= 0 ) CYCLE
                   IF (SIZE(dimlevs) /= SIZE(levs)) CYCLE
                   l_cyc = .FALSE.
                   DO ic = 1, SIZE(dimlevs)
                      IF (dimlevs(ic) /= levs(ic) ) THEN
                         l_cyc=.TRUE.
                         EXIT
                      END IF
                   END DO
                   IF (l_cyc) THEN
                      CYCLE
                   END IF
                   l_mkdim = .FALSE.
                   EXIT
                ENDIF
             END DO
             IF (ASSOCIATED(dimlevs)) THEN
                DEALLOCATE(dimlevs) ; NULLIFY(dimlevs)
             END IF

             IF (l_mkdim) THEN
                curcnt = curcnt + 1
                CALL INT2STR(cntstr,curcnt)
                ! NEW  make every dimension / representation:
                ! REQUIRED to add dimension variables
                CALL NEW_DIMENSION(status, DIMID_VERT &
                     , 'DIM_'//numstr//'LEV'//cntstr, ip, dom_id=np)
                CALL CHANNEL_HALT(substr, status)
                CALL add_dimension_variable(status, DIMID_VERT &
                     , 'DIM_'//numstr//'LEV'//cntstr, levs)
                CALL CHANNEL_HALT(substr, status)
                CALL add_dimension_variable_att(status, DIMID_VERT &
                     , 'DIM_'//numstr//'LEV'//cntstr, 'units', c=TRIM(levunit))
                CALL CHANNEL_HALT(substr, status)
                CALL add_dimension_variable_att(status, DIMID_VERT &
                     , 'DIM_'//numstr//'LEV'//cntstr  &
                     , 'levtype', c=TRIM(dimname))
                CALL CHANNEL_HALT(substr, status)
             END IF
             DEALLOCATE(levs, chlevs)
             ! CHECK ALSO, IF REPRESENTATION EXISTS
             IF (.NOT. l_mkdim) THEN
                CALL get_representation(status &
                     , dim_ids=(/ _RI_XYZ__(DIMID_LON, DIMID_LAT, DIMID_VERT) &
                       , DIMID_UNDEF /) &
                     , reprid=reprid)

                IF ( status /=0 ) THEN
                   l_mkdim = .TRUE.
                   ! extra number required in case DIMID_VERT exists
                   ! was created by another submodel, but 3D representation
                   ! does not exist
                   reprcnt=reprcnt-1
                   CALL INT2STR(cntstr,reprcnt)
                END IF
             END IF
             IF (l_mkdim) THEN
#ifndef ICON
                ! b) make representation
                nseg = gp_nseg

                ALLOCATE(start(nseg,IRANK))
                ALLOCATE(cnt(nseg,IRANK))
                ALLOCATE(meml(nseg,IRANK))
                ALLOCATE(memu(nseg,IRANK))

                start(:,:)=1
                cnt(:,:)  =1
                meml(:,:) =1
                memu(:,:) =1

                CALL NEW_REPRESENTATION(status, reprid &
                     , reprstr//'3D_'//numstr//'LEV'//cntstr &
                     , rank = 3, link = 'xxx-', dctype = DC_GP                 &
                     , dimension_ids = (/ &
                     _RI_XYZ__(DIMID_LON, DIMID_LAT, DIMID_VERT) /)            &
                     , ldimlen       = (/_RI_XYZ__(nproma , ngpblks , AUTO) /) &
                     , output_order  = (/ _IX_XYZ__,_IY_XYZ__, _IZ_XYZ__ /)    &
                     , axis         = repr_def_axes(_RI_XYZ__('X','Y','Z'),'-')&
                     )
                CALL CHANNEL_HALT(substr//'04', status)
                start(:,:) = 1
                cnt(:,:) = 1

                start(:, _IX_XYZ__) = gp_start(:, _IX_XYZN_ )
                cnt(:, _IX_XYZ__)   = gp_cnt(:,_IX_XYZN_)
                meml(:, _IX_XYZ__)  = gp_meml(:,_IX_XYZN_)
                memu(:, _IX_XYZ__)  = gp_memu(:,_IX_XYZN_)

                start(:,_IY_XYZ__) = gp_start(:,_IY_XYZN_)
                cnt(:,_IY_XYZ__)   = gp_cnt(:,_IY_XYZN_)
                meml(:,_IY_XYZ__)  = gp_meml(:,_IY_XYZN_)
                memu(:,_IY_XYZ__)  = gp_memu(:,_IY_XYZN_)

                start(:,_IZ_XYZ__) = gp_start(:,_IZ_XYZN_)
                cnt(:,_IZ_XYZ__)  = ip
                memu(:,_IZ_XYZ__) = ip
                meml(:,_IZ_XYZ__) = gp_meml(:,_IZ_XYZN_)

                CALL SET_REPRESENTATION_DECOMP(status, reprid &
                     , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
                CALL CHANNEL_HALT(substr, status)

                DEALLOCATE(start)
                DEALLOCATE(cnt)
                DEALLOCATE(meml)
                DEALLOCATE(memu)
! def ICON
#else
                CALL NEW_REPRESENTATION(status, reprid                 &
                     , reprstr//'3D_'//numstr//'LEV'//cntstr           &
                     , rank = 3, link = 'xxx-', dctype = DC_GP         &
                     , dimension_ids = (/ DIMID_NCELLS(np), DIMID_VERT &
                     ,                    DIMID_ONE(np) /)             &
                     , ldimlen       = (/ nproma , AUTO                &
                     ,                    p_patch(np)%nblks_c   /)     &
                     , output_order  = (/ 1,3,2 /)                     &
                     , axis          = 'XZY-'                          &
                     , dom_id      = np                                &
                     )
                CALL CHANNEL_HALT(substr//'05', status)
#endif
             END IF
             ! DIAGNOSTIC OUTPUT
             IF (p_parallel_io) THEN
                   WRITE(*,*) '  IMPORT TYPE   : MULTI LEVEL (Nx2D)'
                   WRITE(*,*) '  --> ',ip,' IMPORT LEVELS'
             END IF
             ! ----------------------------------------------------------
          CASE(4)
             DO jv=1, mvar
             ! INDEX REGRIDDING
             !
             ! NUMBER OF LEVELS IN CHANNEL OBJECT
             IF (ALL(IMPORT_SET(i)%IXF(jv,:) == 0)) CYCLE
             ip = (IMPORT_SET(i)%IXF(jv,2) - IMPORT_SET(i)%IXF(jv,1)) + 1
             mlev = ip
             CALL INT2STR(numstr,mlev)

             l_mkdim = .TRUE.
             idxrun = 0
             ALLOCATE(levs(mlev))
             DO ic = 1, mlev
                levs(ic) = IMPORT_SET(i)%IXF(jv,1) + ic - 1
             END DO

             DO
                IF (ASSOCIATED(dimlevs)) THEN
                   DEALLOCATE(dimlevs) ; NULLIFY(dimlevs)
                END IF
                idxrun = idxrun +1
                CALL INT2STR(idxstr, idxrun)
                CALL GET_DIMENSION_INFO(status       &
                     , 'DIM_'//numstr//'IXF'//idxstr &
                     , DIMID_VERT, dlen, dom_id=np)
                IF (status == 905) EXIT ! DIMENSION DOES NOT EXIST
                CALL channel_halt(substr,status)

                IF (dlen == ip) THEN
                   !DIMENSION has same length, test attribute
                   CALL get_dimvar_info(status,DIMID_VERT &
                        , dimlevs, firstname='index')
                   IF (status /= 0 ) CYCLE
                   IF (SIZE(dimlevs) /= SIZE(levs)) CYCLE
                   l_cyc = .FALSE.
                   DO ic = 1, SIZE(dimlevs)
                      IF (dimlevs(ic) /= levs(ic) ) THEN
                         l_cyc=.TRUE.
                         EXIT
                      END IF
                   END DO
                   IF (l_cyc) THEN
                      CYCLE
                   END IF
                   l_mkdim = .FALSE.
                   EXIT
                ENDIF
             END DO
             IF (ASSOCIATED(dimlevs)) THEN
                DEALLOCATE(dimlevs) ; NULLIFY(dimlevs)
             END IF

             IF (l_mkdim) THEN
                curixf = curixf + 1
                CALL INT2STR(cntstr,curixf)
                ! NEW  make every dimension / representation:
                ! REQUIRED to add dimension variables
                CALL NEW_DIMENSION(status, DIMID_VERT &
                     , 'DIM_'//numstr//'IXF'//cntstr, ip, dom_id=np)
                CALL CHANNEL_HALT(substr, status)
                CALL add_dimension_variable(status, DIMID_VERT &
                     , 'index_'//cntstr, levs)
                CALL CHANNEL_HALT(substr, status)
                CALL add_dimension_variable_att(status, DIMID_VERT &
                     , 'index_'//cntstr, 'units', c='1')
                CALL CHANNEL_HALT(substr, status)
             END IF
             DEALLOCATE(levs)
             ! CHECK ALSO, IF REPRESENTATION EXISTS
             IF (.NOT. l_mkdim) THEN
                CALL get_representation(status &
                     , dim_ids=(/ _RI_XYZ__(DIMID_LON, DIMID_LAT, DIMID_VERT) &
                       , DIMID_UNDEF /) &
                     , reprid=reprid)

                IF ( status /=0 ) THEN
                   l_mkdim = .TRUE.
                   ! extra number required in case DIMID_VERT exists
                   ! was created by another submodel, but 3D representation
                   ! does not exist
                   reprcnt=reprcnt-1
                   CALL INT2STR(cntstr,reprcnt)
                END IF
             END IF
             IF (l_mkdim) THEN
#ifndef ICON
                ! b) make representation
                nseg = gp_nseg

                ALLOCATE(start(nseg,IRANK))
                ALLOCATE(cnt(nseg,IRANK))
                ALLOCATE(meml(nseg,IRANK))
                ALLOCATE(memu(nseg,IRANK))

                start(:,:)=1
                cnt(:,:)  =1
                meml(:,:) =1
                memu(:,:) =1

                CALL NEW_REPRESENTATION(status, reprid                       &
                     , reprstr//'3D_'//numstr//'IXF'//cntstr                 &
                     , rank = 3, link = 'xxx-', dctype = DC_GP               &
                     , dimension_ids = (/                                    &
                       _RI_XY_N_(DIMID_LON, DIMID_LAT, DIMID_VERT) /)        &
                     , ldimlen       = (/ _RI_XY_N_(nproma, ngpblks , AUTO) /)&
                     , output_order  = (/ _IX_XY_N_,_IY_XY_N_,_IN_XY_N_ /)   &
                     , axis      = repr_def_axes(_RI_XY_N_('X','Y','N'),'-') &
                     )
                CALL CHANNEL_HALT(substr//'07', status)
                start(:,:) = gp_start(:,:)
                cnt(:,:) = gp_cnt(:,:)
                meml(:,:) = gp_meml(:,:)
                memu(:,:) = gp_memu(:,:)

                start(:,_IX_XY_N_) = gp_start(:,_IX_XYZN_)
                cnt(:,_IX_XY_N_)   = gp_cnt(:,_IX_XYZN_)
                meml(:,_IX_XY_N_)  = gp_meml(:,_IX_XYZN_)
                memu(:,_IX_XY_N_)  = gp_memu(:,_IX_XYZN_)

                start(:,_IY_XY_N_) = gp_start(:,_IY_XYZN_)
                cnt(:,_IY_XY_N_)   = gp_cnt(:,_IY_XYZN_)
                meml(:,_IY_XY_N_)  = gp_meml(:,_IY_XYZN_)
                memu(:,_IY_XY_N_)  = gp_memu(:,_IY_XYZN_)

                start(:,_IN_XY_N_) = 1
                meml(:,_IN_XY_N_) = 1
                cnt(:,_IN_XY_N_)  = mlev
                memu(:,_IN_XY_N_) = mlev

                CALL SET_REPRESENTATION_DECOMP(status, reprid &
                     , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
                CALL CHANNEL_HALT(substr, status)

                DEALLOCATE(start)
                DEALLOCATE(cnt)
                DEALLOCATE(meml)
                DEALLOCATE(memu)
#else
! def ICON
                CALL NEW_REPRESENTATION(status, reprid                 &
                     , reprstr//'3D_'//numstr//'IXF'//cntstr           &
                     , rank = 3, link = 'xxx-', dctype = DC_GP         &
                     , dimension_ids = (/ DIMID_NCELLS(np)             &
                                        , DIMID_ONE(np), DIMID_VERT /) &
                     , ldimlen       = (/ nproma                       &
                                        , p_patch(np)%nblks_c, AUTO /) &
                     , output_order  = (/ 1,2,3 /)                     &
                     , axis          = 'XYN-'                          &
                     , dom_id      = np                                &
                     )
              CALL CHANNEL_HALT(substr//'08', status)
#endif
             END IF
             ! DIAGNOSTIC OUTPUT
             IF (p_parallel_io) THEN
                   WRITE(*,*) '  IMPORT TYPE   : INDEX REGRIDDING'
                   WRITE(*,*) '  --> ',mlev,' INDICES'
             END IF
             ! ----------------------------------------------------------
             END DO
          CASE(5)
             ! 4D field
             !
             ! NUMBER OF PARAMETERS
             CALL INT2STR(numstr,ip)
             CALL GET_REPRESENTATION_INFO(status, reprstr//'4D_'//numstr//'PAR'  &
                  , id=reprid, dom_id=np)
             IF (status == 2003) THEN
                ! in this case the representation was not jet made
                ! A) make vertical dimension
                CALL GET_DIMENSION_INFO(status,'DIM_'//numstr//'PAR' &
                                       , id=DIMID_PAR, dom_id=np)
                IF (status == 905) THEN
                   CALL NEW_DIMENSION(status, DIMID_PAR &
                                    , 'DIM_'//numstr//'PAR', ip, dom_id=np)
                END IF
                CALL CHANNEL_HALT(substr, status)

                ! b) make representation
#ifndef ICON
                nseg = gp_nseg

                ALLOCATE(start(nseg,IRANK))
                ALLOCATE(cnt(nseg,IRANK))
                ALLOCATE(meml(nseg,IRANK))
                ALLOCATE(memu(nseg,IRANK))

                start(:,:)=1
                cnt(:,:)  =1
                meml(:,:) =1
                memu(:,:) =1

                CALL NEW_REPRESENTATION(status, reprid, reprstr//'4D_'//numstr//'PAR'&
                     , rank = 4, link = 'xxxx', dctype = DC_GP                &
                     , dimension_ids = (/ &
                       _RI_XYZN_(DIMID_LON, DIMID_LAT, DIMID_LEV, DIMID_PAR) /)&
                     , ldimlen       = (/&
                        _RI_XYZN_(nproma , ngpblks, AUTO, AUTO)   /)     &
                        , output_order  = (/ _IN_XYZN_, _IX_XYZN_          &
                                           , _IY_XYZN_, _IZ_XYZN_ /)      &
                     , axis          = repr_def_axes('X','Y','Z','N')     &
                     )
                CALL CHANNEL_HALT(substr//'11', status)
                start(:,:) = gp_start(:,:)
                cnt(:,:) = gp_cnt(:,:)

                cnt(:,_IN_XYZN_) = ip
                meml(:,:) = gp_meml(:,:)
                memu(:,:) = gp_memu(:,:)
                memu(:,_IN_XYZN_) = ip

                CALL SET_REPRESENTATION_DECOMP(status, reprid &
                     , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
                CALL CHANNEL_HALT(substr, status)

                DEALLOCATE(start)
                DEALLOCATE(cnt)
                DEALLOCATE(meml)
                DEALLOCATE(memu)
#else
! def ICON
                CALL NEW_REPRESENTATION(status, reprid                       &
                     , reprstr//'4D_'//numstr//'PAR'                            &
                     , rank = 4, link = 'xxxx', dctype = DC_GP               &
                     , dimension_ids = (/ DIMID_NCELLS(np), DIMID_LEV        &
                                        , DIMID_PAR, DIMID_ONE(np) /)        &
                     , ldimlen       = (/ nproma , AUTO, AUTO                &
                                        , p_patch(np)%nblks_c   /)           &
                     , output_order  = (/ 1,4,2,3 /)                         &
                     , axis          = 'XZNY'                                &
                     , dom_id      = np                                    &
                     )
                CALL CHANNEL_HALT(substr//'12', status)
#endif
             ELSE
                CALL CHANNEL_HALT(substr//'13', status)
             ENDIF

             ! DIAGNOSTIC OUTPUT
             IF (p_parallel_io) THEN
                   WRITE(*,*) '  IMPORT TYPE   : 4D REGRIDDING'
                   WRITE(*,*) '  --> ',ip,' IS LENGTH OF PARAMETER DIMENSION'
             END IF
             ! ----------------------------------------------------------
             ! ----------------------------------------------------------
          CASE(6)
             ! 3D field (mlev=1, nparam>1 but no z-Axis)
             !
             ! NUMBER OF PARAMETERS
             CALL INT2STR(numstr,ip)
             CALL GET_REPRESENTATION_INFO(status,reprstr//'3D_'//numstr//'PAR' &
                  , id=reprid, dom_id=np)
             IF (status == 2003) THEN
                ! in this case the representation was not jet made
                ! A) make vertical dimension
                CALL GET_DIMENSION_INFO(status,'DIM_'//numstr//'PAR' &
                                       , id=DIMID_PAR, dom_id=np)
                IF (status == 905) THEN
                   CALL NEW_DIMENSION(status, DIMID_PAR &
                                    , 'DIM_'//numstr//'PAR', ip, dom_id=np)
                END IF
                CALL CHANNEL_HALT(substr, status)

#ifndef ICON
                ! b) make representation
                nseg = gp_nseg

                ALLOCATE(start(nseg,IRANK))
                ALLOCATE(cnt(nseg,IRANK))
                ALLOCATE(meml(nseg,IRANK))
                ALLOCATE(memu(nseg,IRANK))

                start(:,:)=1
                cnt(:,:)  =1
                meml(:,:) =1
                memu(:,:) =1

                CALL NEW_REPRESENTATION(status, reprid                         &
                     , reprstr//'3D_'//numstr//'PAR'                           &
                     , rank = 3, link = 'xxx-', dctype = DC_GP                 &
                     , dimension_ids = (/&
                       _RI_XY_N_(DIMID_LON, DIMID_LAT, DIMID_PAR) /)           &
                     , ldimlen       = (/ _RI_XY_N_(nproma,ngpblks, AUTO)   /) &
                     , output_order  = (/ _IX_XY_N_, _IY_XY_N_, _IN_XY_N_ /)   &
                     , axis        = repr_def_axes(_RI_XY_N_('X','Y','N'),'-') &
                     )
                CALL CHANNEL_HALT(substr//'15', status)
                start(:,:) = gp_start(:,:)
                cnt(:,:) = gp_cnt(:,:)

                cnt(:,_IN_XY_N_) = ip
                meml(:,:) = gp_meml(:,:)
                memu(:,:) = gp_memu(:,:)
                memu(:,_IN_XY_N_) = ip

                CALL SET_REPRESENTATION_DECOMP(status, reprid &
                     , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
                CALL CHANNEL_HALT(substr, status)

                DEALLOCATE(start)
                DEALLOCATE(cnt)
                DEALLOCATE(meml)
                DEALLOCATE(memu)
#else
!ICON
                CALL NEW_REPRESENTATION(status, reprid                       &
                     , reprstr//'3D_'//numstr//'PAR'                         &
                     , rank = 3, link = 'xxx-', dctype = DC_GP               &
                     , dimension_ids = (/ DIMID_NCELLS(np)                   &
                                        , DIMID_ONE(np), DIMID_PAR /)        &
                     , ldimlen       = (/ nproma                             &
                                        , p_patch(np)%nblks_c , AUTO    /)   &
                     , output_order  = (/ 1,2,3 /)                           &
                     , axis          = 'XYN-'                                &
                     , dom_id      = np                                      &
                     )
                CALL CHANNEL_HALT(substr//'16', status)
#endif
             ELSE
                CALL CHANNEL_HALT(substr//'17', status)
             ENDIF

             ! DIAGNOSTIC OUTPUT
             IF (p_parallel_io) THEN
                   WRITE(*,*) '  IMPORT TYPE   : 3D PARAMETER REGRIDDING'
                   WRITE(*,*) '  --> ',ip,' IS LENGTH OF PARAMETER DIMENSION'
             END IF
             ! ----------------------------------------------------------
          CASE DEFAULT
             ! ----------------------------------------------------------
             !
             ! ----------------------------------------------------------
          END SELECT

          variable_loop: DO jv=1, mvar
             !
             IMPORT_SET(i)%nf = IMPORT_SET(i)%nf + 1
             IMPORT_SET(i)%fieldname(jv)%obj = &
                  TRIM(name)//'_'//TRIM(var(jv))
             IF (p_parallel_io) THEN
                WRITE(*,*) '   CHANNEL OBJECT  : '   &
                        ,'(',IMPORT_SET(i)%nf,') '        &
                        ,TRIM(IMPORT_SET(i)%fieldname(jv)%obj)
             END IF
             !
             ! CREATE CHANNEL OBJECTS (2D, 3D, index-regridding, 4D)
             CALL NEW_CHANNEL_OBJECT(status, ch_str &
                  , TRIM(IMPORT_SET(i)%fieldname(jv)%obj) &
                  , p4=IMPORT_SET(i)%field(jv)%ptr &
                  , reprid=reprid                  &
                  , lstatic=lstatic                &
                  , lrestreq = lrestart            )
             CALL CHANNEL_HALT(substr, status)

             IF (ASSOCIATED(att)) THEN
                IF (ASSOCIATED(att(jv)%ptr)) THEN
                   CALL ADD_VAR_ATTRIBUTES(status, ch_str &
                        , TRIM(IMPORT_SET(i)%fieldname(jv)%obj) &
                        , att(jv)%ptr )
                   CALL CHANNEL_HALT(substr, status)
                   DEALLOCATE(att(jv)%ptr) ; NULLIFY(att(jv)%ptr)
                END IF
             END IF
             CALL NEW_ATTRIBUTE(status, ch_str &
                  ,  TRIM(IMPORT_SET(i)%fieldname(jv)%obj) &
                  , 'mmig_event_name', c=TRIM(name) )
             !
          END DO variable_loop

          IF (ASSOCIATED(att)) THEN
             DEALLOCATE(att) ; NULLIFY(att)
          END IF

       ELSE ! NOT OK
          IF (p_parallel_io) THEN
             WRITE(*,*) 'ERROR:'
             WRITE(*,*) '   REGRIDDING OF RGT-EVENT ',i,' FAILED !'
             CALL ERROR_BI(' ',substr)
          END IF
       END IF      ! IF (lok)
       !

       ! CLEAN UP
       IF (ASSOCIATED(dat)) DEALLOCATE(dat)
       NULLIFY(dat)
       IF (ASSOCIATED(var)) DEALLOCATE(var)
       NULLIFY(var)

       IF (p_parallel_io) THEN
          WRITE(*,*) '========================================================'
       END IF

    END DO import_set_loop

    CALL END_MESSAGE_BI(ch_str, 'MEMORY SETUP',substr)

  END SUBROUTINE IMPORT_GRID_INIT_MEMORY
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE IMPORT_GRID_INIT_TRACER

    IMPLICIT NONE

    CALL IMPORT_GRID_REGRID(L_INIT=.TRUE.)

  END SUBROUTINE IMPORT_GRID_INIT_TRACER
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE IMPORT_GRID_GLOBAL_START

    IMPLICIT NONE

    CALL IMPORT_GRID_REGRID(L_INIT=.FALSE.)

  END SUBROUTINE IMPORT_GRID_GLOBAL_START
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE IMPORT_GRID_REGRID(L_INIT)

    USE messy_main_mpi_bi,               ONLY: p_parallel_io
#if defined(COSMO) || defined(MESSYDWARF)
    USE messy_main_data_bi,              ONLY: basemodstr => modstr
#endif
    USE messy_main_grid_def_mem_bi,      ONLY: BASEGRID_ID
    USE messy_main_channel_error_bi,     ONLY: CHANNEL_HALT
    USE messy_main_import_grid_tools_bi, ONLY: RGTOOL_BI_READ_NCVAR  &
                                             , RGTOOL_BI_READ_NCFILE &
                                             , RGTEVENT_STATUS &
                                             , RGTEVENT_FROM_CHANNEL
    USE messy_main_grid_netcdf,          ONLY: GRD_MAXSTRLEN &
                                             , QDEF_NCVAR
    USE messy_main_channel_mem,          ONLY: dom_curid
    USE messy_main_channel,              ONLY: GET_CHANNEL_OBJECT          &
                                             , t_chaobj_cpl, STRLEN_OBJECT &
                                             , GET_CHANNEL_OBJECT_INFO
    USE messy_main_channel_repr,         ONLY: GET_REPRESENTATION_INFO
    USE messy_main_tools,                ONLY: lcase
    USE messy_main_timer,                ONLY: lresume
    USE messy_main_grid,                 ONLY: LOCATE_GEOHYBGRID           &
                                             , INIT_GEOHYBGRID             &
                                             , GEOHYBGRID_UPDATE_PRESS
    USE messy_main_grid_tools,           ONLY: RGTOOL_CONVERT_DAT2PREDEFVAR

    IMPLICIT NONE

    ! I/O
    LOGICAL, INTENT(IN) :: L_INIT

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'import_grid_regrid'
    INTEGER :: i, jv, iz
    LOGICAL :: l_event
    INTEGER :: tnc
    REAL(dp),                     DIMENSION(:,:,:,:,:), POINTER :: dat => NULL()
    CHARACTER(LEN=GRD_MAXSTRLEN), DIMENSION(:),         POINTER :: var => NULL()
    INTEGER :: mvar
    LOGICAL :: lok
    INTEGER                     :: np
    CHARACTER(LEN=4)            :: caxis   = ' '
    TYPE(t_geohybgrid)          :: xgrid
    INTEGER                     :: status
    INTEGER                     :: reprid = 0
    TYPE(t_chaobj_cpl)          :: pressiobj
    TYPE(t_chaobj_cpl)          :: pressmobj
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: press4d => NULL()
    LOGICAL,  DIMENSION(:),   ALLOCATABLE :: lassign
    LOGICAL                     :: lread  = .FALSE.
    LOGICAL                     :: lrread = .FALSE.

    INTEGER                       :: dom_idx
    CHARACTER(LEN=STRLEN_OBJECT)  :: name
    INTEGER                       :: LOOP_START, LOOP_END

#ifdef COSMO
    pressiobj%CHA = basemodstr
    pressiobj%OBJ = 'pressi'
    pressmobj%CHA = basemodstr
    pressmobj%OBJ = 'press'
#elif defined(ECHAM5) || defined(CESM1)
    ! NOT REQUIRED, pressure grid
#elif defined(ICON)
    ! variables added to varlist in atm_dyn_iconam/mo_nonhydro_state.f90
    pressmobj%CHA = 'nh_state_diag'
    pressmobj%OBJ = 'pres'
    pressiobj%CHA = 'nh_state_diag'
    pressiobj%OBJ = 'pres_ifc'
#endif

    ! set patch
    np  = dom_curid

    IF (NIMPALL_SET == 0) RETURN

    IF (L_INIT) THEN
       IF (NIMPINIT_SET == 0) RETURN
       LOOP_START = NIMPORT_SET + 1
       LOOP_END   = NIMPALL_SET
    ELSE
       IF (NIMPORT_SET == 0)  RETURN
       LOOP_START = 1
       LOOP_END   = NIMPORT_SET

       ! SET INTERNAL COUNTER FROM CHANNEL OBJECT
       IF (lresume) CALL RGTEVENT_FROM_CHANNEL(rgt)
    END IF

    IF (L_req3Dipol(np)) THEN
       ! UPDATE PRESSURE FIELD OF BASEGRID IF REQUIRED
       CALL LOCATE_GEOHYBGRID(status, BASEGRID_ID(np), grid=xgrid)
       IF (status /= 0) CALL error_bi('BASEGRID NOT LOCATABLE',substr)

       ! GET PRESSURE FIELD, if base grid is height grid
       ! => ECHAM and CESM do not need these fields (xgrid%lpres = F)
       ! => COSMO and ICON require it
       IF (QDEF_NCVAR(xgrid%pressm)) THEN ! is pressm defined?
          ! GET PRESSURE FIELD
          CALL GET_CHANNEL_OBJECT(status, TRIM(pressmobj%cha) &
               , TRIM(pressmobj%OBJ) &
               , p4=press4d, dom_id=np)
          CALL CHANNEL_HALT(substr, status)
          CALL GET_CHANNEL_OBJECT_INFO(status &
               ,TRIM(pressmobj%cha),TRIM(pressmobj%obj), reprid=reprid &
               ,dom_id=np )
          CALL CHANNEL_HALT(substr, status)
          CALL GET_REPRESENTATION_INFO(status, ' ', ID=reprid, axis=caxis)
          CALL CHANNEL_HALT(substr, status)
          CALL lcase(caxis)
          CALL RGTOOL_CONVERT_DAT2PREDEFVAR(status, xgrid%pressm &
               , press4d, caxis, 'xyzn')

          !NULLIFY(press4d)
          CALL GEOHYBGRID_UPDATE_PRESS(status, BASEGRID_ID(np) &
               , pressm=xgrid%pressm )
          NULLIFY(press4d)
       END IF

       IF (QDEF_NCVAR(xgrid%pressi)) THEN ! is pressi defined?
          ! GET PRESSURE FIELD
          CALL GET_CHANNEL_OBJECT(status, TRIM(pressiobj%cha) &
               , TRIM(pressiobj%OBJ) &
               , p4=press4d, dom_id=np)
          CALL CHANNEL_HALT(substr, status)
          CALL GET_CHANNEL_OBJECT_INFO(status &
               , TRIM(pressiobj%cha),TRIM(pressiobj%obj), reprid=reprid &
               , dom_id=np )
          CALL CHANNEL_HALT(substr, status)
          CALL GET_REPRESENTATION_INFO(status, ' ', ID=reprid, axis=caxis)
          CALL CHANNEL_HALT(substr, status)
          CALL lcase(caxis)
          CALL RGTOOL_CONVERT_DAT2PREDEFVAR(status, xgrid%pressi &
               , press4d, caxis, 'xyzn')
          CALL GEOHYBGRID_UPDATE_PRESS(status, BASEGRID_ID(np) &
               , pressi=xgrid%pressi )
          NULLIFY(press4d)
       END IF
       CALL INIT_GEOHYBGRID(xgrid)
    ENDIF

    rgt_set_loop: DO i=LOOP_START, LOOP_END

       IF (L_INIT) THEN
          dom_idx = RGTINIT(i-NIMPORT_SET)%domain_idx
          name    = RGTINIT(i-NIMPORT_SET)%name
          IF (ALLOCATED(lassign)) DEALLOCATE(lassign)
          ALLOCATE(lassign(SIZE(IMPORT_SET(i)%fieldname)))
          lassign(:) = .TRUE.
       ELSE
          dom_idx = RGT(i)%domain_idx
          name    = RGT(i)%io%cnt%name
       END IF

       IF (.NOT. IMPORT_SET(i)%lok) CYCLE     ! IMPORT SET NOT OK
       IF (dom_idx /= dom_curid) CYCLE

       lread = .FALSE.
       IF (lresume) THEN
          DO jv = 1, SIZE(IMPORT_SET(i)%fieldname)
             lrread = .FALSE.
             CALL GET_CHANNEL_OBJECT_INFO(status          &
                  , TRIM(IMPORT_SET(i)%fieldname(jv)%cha) &
                  , TRIM(IMPORT_SET(i)%fieldname(jv)%obj) &
                  , lrestart_read=lrread)
             lread = lread .OR. (.NOT. lrread)
             IF (L_INIT) lassign(jv) = .NOT. lrread
          END DO
       END IF

       IF (.NOT. L_INIT) THEN
          CALL RGTEVENT_STATUS(l_event, tnc, RGT, ix = i &
               , lstop=.true., lftrig=lread)

          IF (.NOT. l_event) CYCLE             ! EVENT NOT TRIGGERED
       ELSE
          IF (.NOT. ANY(lassign)) THEN
             DEALLOCATE(lassign)
             CYCLE
          END IF
          tnc = RGTINIT(i-NIMPORT_SET)%tstep
       END IF

       IF (p_parallel_io) THEN
          WRITE(*,*) '========================================================'
          WRITE(*,*) 'UPDATING IMPORT-SET ',i,' ...'
       END IF

       ! define output grid: DEFAULT: base model grid
       ! TODO: add variable output grid

       SELECT CASE(IMPORT_SET(i)%RGREAD)
       CASE(RGREAD_NCFILE)
          !
          CALL RGTOOL_BI_READ_NCFILE( TRIM(IMPORT_SET(i)%nml)   &
               , TRIM(IMPORT_SET(i)%file), tnc, dat             &
               , IMPORT_SET(i)%GTRF_TYPE                        &
               , IMPORT_SET(i)%IXF                              &
               , IMPORT_SET(i)%ogrid_id, IMPORT_SET(i)%igrid_id &
               , var, SDID=IMPORT_SET(i)%SCRIP_DATA_ID          &
               , lrg = .true., lrgz=( (IMPORT_SET(i)%itype==2)  &
                                 .OR. (IMPORT_SET(i)%itype==5)) &
               , lok = lok                                      &
               , mnc=IMPORT_SET(i)%mnc )
          !
       CASE(RGREAD_NCVAR)
          !
          CALL RGTOOL_BI_READ_NCVAR( TRIM(IMPORT_SET(i)%nml)    &
               , TRIM(IMPORT_SET(i)%var), tnc, dat              &
               , IMPORT_SET(i)%GTRF_TYPE                        &
               , IMPORT_SET(i)%IXF                              &
               , IMPORT_SET(i)%ogrid_id, IMPORT_SET(i)%igrid_id &
               , SDID=IMPORT_SET(i)%SCRIP_DATA_ID               &
               , lrg = .true., lrgz=((IMPORT_SET(i)%itype==2)   &
                                .OR. (IMPORT_SET(i)%itype==5))  &
               , lok = lok                                      &
               , mnc=IMPORT_SET(i)%mnc )
          ALLOCATE(var(1))
          var(1) = TRIM(IMPORT_SET(i)%var)
          !
       CASE DEFAULT
          !
          ! SHOULD BE NEVER REACHED
          CALL ERROR_BI('UNKNOWN RGREAD METHOD',substr)
          !
       END SELECT

       IF (.NOT. lok) CALL ERROR_BI('REGRIDDING NOT SUCCESSFUL !',substr)

       mvar = SIZE(var)
       IF (mvar /= IMPORT_SET(i)%nvar) &
            CALL ERROR_BI('NUMBER OF VARIABLES CHANGED !',substr)

       IF (p_parallel_io) THEN
          WRITE(*,*) ' RGTEVENT:   ',TRIM(name)
       END IF

       variable_loop: DO jv=1, mvar

          IF (L_INIT) THEN
             IF (.NOT. lassign(jv)) CYCLE
          END IF

          ! NO CHANNEL OBJECT PRESENT
          IF (.NOT. ASSOCIATED(IMPORT_SET(i)%field(jv)%ptr)) CYCLE

          IF (p_parallel_io) THEN
             WRITE(*,*) '   VARIABLE: ',TRIM(var(jv))
          END IF

          SELECT CASE(IMPORT_SET(i)%itype)
          CASE(0) ! RAW / no decomposition
             IF (p_parallel_io) THEN
                WRITE(*,*) '   ... RAW'
             END IF
             IF (SIZE(dat,3) > 1) THEN
                IMPORT_SET(i)%field(jv)%ptr(:,:,:,1) = REAL(dat(:,:,:,1,jv),dp)
             ELSE
                IMPORT_SET(i)%field(jv)%ptr(:,:,:,1) = REAL(dat(:,:,1,:,jv),dp)
             ENDIF
          CASE(1) ! 2D
             IF (p_parallel_io) THEN
                WRITE(*,*) '   ... 2D'
             END IF
             IF (SIZE(dat,_IZ_XYZN_) /= 1 .OR. SIZE(dat,_IN_XYZN_) /= 1) THEN
                CALL ERROR_BI( &
                     'PARAMETER/LEVEL DIMENSION OF 2D FIELD IS > 1' &
                     , substr)
             END IF
             IMPORT_SET(i)%field(jv)%ptr(_RI_XYZ__(:,:,:),1) = &
                  REAL(dat(_RI_XYZN_(:,:,:,1),jv),dp)

          CASE(2) ! 3D
             IF (p_parallel_io) THEN
                WRITE(*,*) '   ... 3D'
             END IF
             IMPORT_SET(i)%field(jv)%ptr(:,:,:,1) = &
                  REAL(dat(_RI_XYZN_(:,:,:,1),jv),dp)

          CASE(3) ! Nx2D

             ! (jv = mvar = 1; since NCVAR processing !)
             IF (p_parallel_io) THEN
                WRITE(*,*) '   ... Nx2D (N=',SIZE(dat,_IN_XYZN_),')'
             END IF
             DO iz = 1, SIZE(dat,_IN_XYZN_)
                IMPORT_SET(i)%field(jv)%ptr(_RI_XYZ__(:,:,iz),1)  = &
                   REAL(dat(_RI_XYZN_(:,:,1,iz),jv),dp)
             END DO
          CASE(4) ! INDEX REGRIDDING
             IF (p_parallel_io) THEN
                WRITE(*,*) '   ... IXF (N=',SIZE(dat,_IN_XYZN_),')'
             END IF
             IMPORT_SET(i)%field(jv)%ptr(:,:,:,1) = 0.0_dp
             IMPORT_SET(i)%field(jv)%ptr(&
                  _RI_XY_N_(:,:,1:SIZE(dat,_IN_XYZN_)),1) = &
                  dat(_RI_XYZN_(:,:,1,:),jv)
          CASE(5) ! 4D
             IF (p_parallel_io) THEN
                WRITE(*,*) '   ... 4D (N=',SIZE(dat,_IN_XYZN_),')'
             END IF
             IMPORT_SET(i)%field(jv)%ptr(:,:,:,:) = dat(:,:,:,:,jv)
          CASE(6) ! 3D x-y-npara
             IF (p_parallel_io) THEN
                WRITE(*,*) '   ... 3D (N=',SIZE(dat,_IN_XYZN_),')'
             END IF
             IMPORT_SET(i)%field(jv)%ptr(&
                  _RI_XY_N_(:,:,1:SIZE(dat,_IN_XYZN_)),1) = &
                  dat(_RI_XYZN_(:,:,1,:),jv)
          CASE DEFAULT

          END SELECT

       END DO variable_loop

       ! CLEAN UP
       IF (ASSOCIATED(dat)) DEALLOCATE(dat)
       NULLIFY(dat)
       IF (ASSOCIATED(var)) DEALLOCATE(var)
       NULLIFY(var)
       IF (ALLOCATED(lassign)) DEALLOCATE(lassign)

!       IF (p_parallel_io) THEN
          WRITE(*,*) 'END UPDATING IMPORT-SET ',i,' !'
          WRITE(*,*) '========================================================'
!       END IF

    END DO rgt_set_loop

  END SUBROUTINE IMPORT_GRID_REGRID
! ------------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE IMPORT_GRID_FREE_MEMORY

    USE messy_main_grid_netcdf,          ONLY: init_multinc

    IMPLICIT NONE

    ! LOCAL
    INTEGER :: i
    CHARACTER(LEN=*),PARAMETER :: substr = 'IMPORT_GRID_FREE_MEMORY'

    CALL START_MESSAGE_BI(submodstr, 'GLOBAL SETUP',substr)

    DO i=1, NIMPALL_SET

       ! FIELD
       IF (ASSOCIATED(IMPORT_SET(i)%field)) &
               DEALLOCATE(IMPORT_SET(i)%field)
       ! Z
       IF (ASSOCIATED(IMPORT_SET(i)%z)) &
            DEALLOCATE(IMPORT_SET(i)%z)
       ! P
       IF (ASSOCIATED(IMPORT_SET(i)%p)) &
            DEALLOCATE(IMPORT_SET(i)%p)

       ! MNC
       CALL init_multinc(IMPORT_SET(i)%mnc)

       IF (ASSOCIATED(IMPORT_SET(i)%IXF)) &
            DEALLOCATE(IMPORT_SET(i)%IXF)
    END DO

    DEALLOCATE(L_req3Dipol)

    CALL END_MESSAGE_BI(submodstr, 'GLOBAL SETUP',substr)

  END SUBROUTINE IMPORT_GRID_FREE_MEMORY
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE ADD_VAR_ATTRIBUTES(status, cha, obj, att)

    USE messy_main_channel_attributes, ONLY: t_attribute &
                                           , TYPE_INTEGER, TYPE_STRING &
                                           , TYPE_REAL_DP
    USE messy_main_channel,            ONLY: NEW_ATTRIBUTE

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER,           INTENT(OUT) :: status
    CHARACTER(LEN=*),  INTENT(IN)  :: cha
    CHARACTER(LEN=*),  INTENT(IN)  :: obj
    ! INTENT(IN)
    TYPE(t_attribute), DIMENSION(:), POINTER     :: att ! attribute list

    ! LOCAL
    INTEGER :: n, i

    status = 0

    IF (.NOT. ASSOCIATED(att)) RETURN

    n = SIZE(att)
    DO i=1, n

       SELECT CASE(att(i)%type)
       CASE(TYPE_INTEGER)
          CALL NEW_ATTRIBUTE(status, TRIM(cha), TRIM(obj) &
               , TRIM(att(i)%name), i=att(i)%i)
       CASE(TYPE_STRING)
          CALL NEW_ATTRIBUTE(status, TRIM(cha), TRIM(obj) &
               , TRIM(att(i)%name), c=att(i)%c)
       CASE(TYPE_REAL_DP)
          CALL NEW_ATTRIBUTE(status, TRIM(cha), TRIM(obj) &
               , TRIM(att(i)%name), r=att(i)%r)
       END SELECT
       IF (status /= 0) RETURN

    END DO

    status = 0

  END SUBROUTINE ADD_VAR_ATTRIBUTES
  ! ----------------------------------------------------------------------

! BLANK
#endif

#if defined(MBM_QBO)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: import_grid_initialize
  PUBLIC :: import_grid_init_memory
  PUBLIC :: import_grid_init_tracer
  PUBLIC :: import_grid_global_start
  PUBLIC :: import_grid_free_memory

CONTAINS

  SUBROUTINE IMPORT_GRID_INITIALIZE

    IMPLICIT NONE

  END SUBROUTINE IMPORT_GRID_INITIALIZE

  SUBROUTINE IMPORT_GRID_INIT_MEMORY

    IMPLICIT NONE

  END SUBROUTINE IMPORT_GRID_INIT_MEMORY

  SUBROUTINE IMPORT_GRID_INIT_TRACER

    IMPLICIT NONE

  END SUBROUTINE IMPORT_GRID_INIT_TRACER

  SUBROUTINE IMPORT_GRID_GLOBAL_START

    IMPLICIT NONE

  END SUBROUTINE IMPORT_GRID_GLOBAL_START

  SUBROUTINE IMPORT_GRID_FREE_MEMORY

    IMPLICIT NONE

  END SUBROUTINE IMPORT_GRID_FREE_MEMORY

#endif

!*************************************************************************
END MODULE MESSY_MAIN_IMPORT_GRID_BI
!*************************************************************************
