!******************************************************************************
MODULE messy_main_grid_bi
!******************************************************************************

  USE messy_main_grid,            ONLY: t_geohybgrid
  USE messy_main_grid_def_mem_bi, ONLY: BASEGRID_ID
  USE messy_main_blather_bi,      ONLY: start_message_bi, end_message_bi
  USE messy_main_constants_mem,   ONLY: dp

  IMPLICIT NONE
  PRIVATE
  SAVE

  ! DEFINE MODEL BASEGRID
  ! DUE to cyclic dependency in COSMO/MESSy
  !     tracer_bi => grid_bi => src_input => src_tracer => tracer_bi
  ! BASEGRID_ID can not be defined here:
  ! There are two solutions possible:
  !   a) create grid_mem_bi  (only for this integer
  !   b) move BASEGRID_ID TO DATA_BI
  ! b) is enacted here
  !INTEGER,                     PUBLIC :: BASEGRID_ID

  ! SUBROUTINES
  PUBLIC :: main_grid_setup
  PUBLIC :: main_grid_initialize
  PUBLIC :: main_grid_init_memory
#if defined(COSMO)
  PUBLIC :: main_grid_init_coupling
#endif
  PUBLIC :: main_grid_read_restart
  PUBLIC :: main_grid_global_start
  PUBLIC :: main_grid_local_start
#if defined(ECHAM5)
  PUBLIC :: main_grid_radiation
#endif
  PUBLIC :: main_grid_vdiff
  PUBLIC :: main_grid_free_memory
  PUBLIC :: p_bcast_grid
#ifdef ICON
  PUBLIC :: main_grid_set_domain
  PUBLIC :: main_grid_set_domain_local
#endif

CONTAINS

  ! ---------------------------------------------------------------------------
  SUBROUTINE p_bcast_grid(grid, proc)

    USE messy_main_mpi_bi,         ONLY: P_BCAST, p_pe
    USE messy_main_grid_netcdf_bi, ONLY: P_BCAST_NCVAR, P_BCAST_NCATT
    USE messy_main_grid,           ONLY: INIT_GEOHYBGRID

    IMPLICIT NONE

    ! I/O
    TYPE(t_geohybgrid), INTENT(INOUT) :: grid
    INTEGER           , INTENT(IN)    :: proc
    ! LOCAL
    INTEGER                           :: i

    IF (p_pe /= proc) THEN
       CALL INIT_GEOHYBGRID(grid)
    ENDIF

    CALL P_BCAST(grid%name, proc)
    CALL P_BCAST(grid%ID,   proc)
    CALL P_BCAST(grid%file, proc)
    CALL P_BCAST(grid%t,    proc)

    CALL P_BCAST_NCATT(grid%att, proc)

    DO i = 1, 4
       CALL P_BCAST(grid%ranges(i,1), proc)
       CALL P_BCAST(grid%ranges(i,2), proc)
    END DO

    CALL P_BCAST(grid%corners,           proc)
    CALL P_BCAST(grid%minmaxlonlat(:,:), proc)
    CALL P_BCAST(grid%start(:),          proc)
    CALL P_BCAST(grid%count(:),          proc)
    CALL P_BCAST(grid%lonc,              proc)

    CALL P_BCAST_NCVAR(grid%lonm,  proc)
    CALL P_BCAST_NCVAR(grid%latm,  proc)
    CALL P_BCAST_NCVAR(grid%hyam,  proc)
    CALL P_BCAST_NCVAR(grid%hybm,  proc)
    CALL P_BCAST_NCVAR(grid%timem, proc)

    CALL P_BCAST_NCVAR(grid%loni,  proc)
    CALL P_BCAST_NCVAR(grid%lati,  proc)
    CALL P_BCAST_NCVAR(grid%hyai,  proc)
    CALL P_BCAST_NCVAR(grid%hybi,  proc)
    CALL P_BCAST_NCVAR(grid%timei, proc)

    CALL P_BCAST_NCVAR(grid%ps,    proc)
    CALL P_BCAST_NCVAR(grid%p0,    proc)

    CALL P_BCAST(grid%clonc,       proc)

    CALL P_BCAST_NCVAR(grid%clonm,  proc)
    CALL P_BCAST_NCVAR(grid%clatm,  proc)
    CALL P_BCAST_NCVAR(grid%cloni,  proc)
    CALL P_BCAST_NCVAR(grid%clati,  proc)

    CALL P_BCAST(grid%rlonc,       proc)

    CALL P_BCAST_NCVAR(grid%rlonm,  proc)
    CALL P_BCAST_NCVAR(grid%rlatm,  proc)
    CALL P_BCAST_NCVAR(grid%rloni,  proc)
    CALL P_BCAST_NCVAR(grid%rlati,  proc)

    CALL P_BCAST_NCVAR(grid%pollon, proc)
    CALL P_BCAST_NCVAR(grid%pollat, proc)
    CALL P_BCAST_NCVAR(grid%polgam, proc)

    CALL P_BCAST_NCVAR(grid%ulonm,  proc)
    CALL P_BCAST_NCVAR(grid%ulatm,  proc)
    CALL P_BCAST_NCVAR(grid%uloni,  proc)
    CALL P_BCAST_NCVAR(grid%ulati,  proc)

    CALL P_BCAST_NCVAR(grid%ps,     proc)
    CALL P_BCAST_NCVAR(grid%p0,     proc)

    CALL P_BCAST_NCVAR(grid%imask,  proc)

    CALL P_BCAST_NCVAR(grid%pressm, proc)
    CALL P_BCAST_NCVAR(grid%pressi, proc)

  END SUBROUTINE p_bcast_grid
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE main_grid_setup

#if defined(BLANK)
    USE messy_main_grid_def_bi, ONLY: main_grid_def_setup

    IMPLICIT NONE

    CALL main_grid_def_setup
#endif

  END SUBROUTINE main_grid_setup
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE main_grid_initialize

    USE messy_main_grid_def_bi, ONLY: main_grid_def_initialize

    IMPLICIT NONE

    CALL main_grid_def_initialize

  END SUBROUTINE main_grid_initialize
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE main_grid_init_memory(flag)


    USE messy_main_grid_def_bi, ONLY: main_grid_def_init_memory
#ifdef COSMO
    USE messy_main_grid_def_bi, ONLY: GRID_DEF_READ_VERTAXIS
#endif
    USE messy_main_blather_bi,  ONLY: ERROR_BI
    USE messy_main_channel_bi,  ONLY: n_dom
    USE messy_main_grid,        ONLY: NEW_GEOHYBGRID, GRID_ERROR &
                                    , INIT_GEOHYBGRID
    USE messy_main_channel_mem, ONLY: dom_curid

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)         :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_grid_init_memory'
    TYPE(t_geohybgrid)          :: bgrid     ! base model grid
    INTEGER                     :: status
    LOGICAL                     :: lok
    INTEGER                     :: np

    np = dom_curid

#ifndef MBM_MPIOM
    IF (.NOT. ALLOCATED(BASEGRID_ID)) THEN
       ALLOCATE(BASEGRID_ID(0:n_dom))
       BASEGRID_ID(:) = -99
    END IF

    SELECT CASE(flag)
    CASE (1)
#ifdef COSMO
       CALL GRID_DEF_READ_VERTAXIS
#endif
       CALL main_grid_def_init_memory

    CASE(2)
       CALL start_message_bi('DEFINING BASEMODEL GRID', 'INIT_MEMORY',substr)
       !  define base model grid:
       CALL DEFINE_BASEMODEL_GEOHYBGRID(bgrid, lok)
       IF (.NOT. lok) &
            CALL ERROR_BI ('definition of basemodel grid failed', substr)

       ! copy grid into list of grids
       CALL NEW_GEOHYBGRID(status, BASEGRID_ID(np), bgrid)
       IF (status /= 0 .AND. status /= 01) &
            CALL ERROR_BI(grid_error(status), substr)

       ! free memory
       CALL INIT_GEOHYBGRID(bgrid)

       CALL end_message_bi('DEFINING BASEMODEL GRID', 'INIT_MEMORY',substr)

    END SELECT

  CONTAINS

#endif

#ifdef ECHAM5
    ! -------------------------------------------------------------------------
    SUBROUTINE DEFINE_BASEMODEL_GEOHYBGRID(g, ok)

      USE messy_main_grid_def_mem_bi, ONLY: nvclev, vct, nlon, ngl, nhgl, nlev
      USE messy_main_grid_def_bi,     ONLY: gl_gmu
      USE messy_main_data_bi,         ONLY:  aps

      USE messy_main_timer,         ONLY: YEAR_START, MONTH_START, DAY_START &
                                        , HOUR_START, MINUTE_START           &
                                        , SECOND_START                       &
                                        , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
      USE messy_main_constants_mem, ONLY: api=>pi, sp
      USE messy_main_timer,         ONLY: TIME_SPAN_D

      USE messy_main_grid_trafo,   ONLY: RGEMPTY, COMPLETE_GEOHYBGRID
      USE messy_main_grid_netcdf,  ONLY: ERRMSG  &
                                       , NULL_DIMID, NULL_VARID &
                                       , VTYPE_REAL, VTYPE_DOUBLE &
                                       , nf90_float  &
                                       , POSITION, INIT_NARRAY, ADD_NCATT
      USE messy_main_grid,         ONLY: t_geohybgrid, INIT_GEOHYBGRID

      IMPLICIT NONE

      ! I/O
      TYPE (t_geohybgrid), INTENT(OUT) :: g
      LOGICAL          ,   INTENT(OUT) :: ok

      ! LOCAL
      REAL(DP)                 :: dts
      INTEGER                  :: i,j
      INTEGER                  :: status
      CHARACTER(LEN=100)       :: tunit

      ! INIT
      CALL INIT_GEOHYBGRID(g)

      g%name = 'ECHAM5'
      g%file = ''                              ! Filename
      g%t    = 0                               ! time step
      g%corners = 4                            ! number of corners

      ! LONGITUDE (MID) ...
      g%lonm%name  = 'lon'
      g%lonm%id    = NULL_VARID
      g%lonm%xtype = NF90_FLOAT
      ! ... dimensions
      g%lonm%ndims = 1
      ALLOCATE(g%lonm%dim(g%lonm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%lonm%dim(1)%name  = 'lon'
      g%lonm%dim(1)%id    = NULL_DIMID
      g%lonm%dim(1)%len   = nlon
      g%lonm%dim(1)%fuid  = .false.
      g%lonm%dim(1)%varid = NULL_VARID
      g%lonc = .TRUE. ! Axis is circular
      ! ... data
      CALL INIT_NARRAY(g%lonm%dat, g%lonm%ndims, (/g%lonm%dim(1)%len/) &
           ,VTYPE_REAL)
      DO i=1, nlon
         g%lonm%dat%vr(i) = 360. * REAL(i-1) / REAL(nlon)
      END DO
      ! ... attributes
      CALL ADD_NCATT(g%lonm, 'long_name', vs='longitude')
      CALL ADD_NCATT(g%lonm, 'units', vs='degrees_east')
      g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
      g%ranges(1,2) = RGEMPTY

      ! LATITUDE (MID) ...
      g%latm%name  = 'lat'
      g%latm%id    = NULL_VARID
      g%latm%xtype = NF90_FLOAT
      ! ... dimensions
      g%latm%ndims = 1
      ALLOCATE(g%latm%dim(g%latm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
      g%latm%dim(1)%name  = 'lat'
      g%latm%dim(1)%id    = NULL_DIMID
      g%latm%dim(1)%len   = ngl
      g%latm%dim(1)%fuid  = .false.
      g%latm%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%latm%dat, g%latm%ndims, (/g%latm%dim(1)%len/) &
           ,VTYPE_REAL)

      DO j=1,nhgl
         g%latm%dat%vr(j)=REAL(ASIN(gl_gmu(j))*180./api,sp)
         g%latm%dat%vr(ngl+1-j)=-g%latm%dat%vr(j)
      END DO

      ! ... attributes
      CALL ADD_NCATT(g%latm, 'long_name', vs='latitude')
      CALL ADD_NCATT(g%latm, 'units', vs='degrees_north')
      g%ranges(2,1) = -90.0_dp
      g%ranges(2,2) =  90.0_dp

      ! HYBRID-A-COEFFICIENTS (INTERFACES) ...
      g%hyai%name  = 'hyai'
      g%hyai%id    = NULL_VARID
      g%hyai%xtype = NF90_FLOAT
      ! ... dimensions
      g%hyai%ndims = 1
      ALLOCATE(g%hyai%dim(g%hyai%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
      g%hyai%dim(1)%name  = 'ilev'
      g%hyai%dim(1)%id    = NULL_DIMID
      g%hyai%dim(1)%len   = nlev+1
      g%hyai%dim(1)%fuid  = .false.
      g%hyai%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%hyai%dat, g%hyai%ndims, (/g%hyai%dim(1)%len/) &
           ,VTYPE_REAL)
      g%hyai%dat%vr(:) = REAL(vct(1:nvclev),sp)
      ! ... attributes
      CALL ADD_NCATT(g%hyai, 'long_name'                              &
           ,vs='hybrid-A-coefficients at layer interfaces')
      CALL ADD_NCATT(g%hyai, 'units', vs='1')
      g%ranges(3,1) = 0.0_dp
      g%ranges(3,2) = 0.0_dp

      ! HYBRID-B-COEFFICIENTS (INTERFACES) ...
      g%hybi%name  = 'hybi'
      g%hybi%id    = NULL_VARID
      g%hybi%xtype = NF90_FLOAT
      ! ... dimensions
      g%hybi%ndims = 1
      ALLOCATE(g%hybi%dim(g%hybi%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,4)
      g%hybi%dim(1) = g%hyai%dim(1)
      ! ... data
      CALL INIT_NARRAY(g%hybi%dat, g%hybi%ndims, (/g%hybi%dim(1)%len/) &
           ,VTYPE_REAL)
      g%hybi%dat%vr(:) = REAL(vct(nvclev+1:2*nvclev),sp)
      ! ... attributes
      CALL ADD_NCATT(g%hybi, 'long_name'                              &
           ,vs='hybrid-B-coefficients at layer interfaces')
      CALL ADD_NCATT(g%hybi, 'units', vs='1')
      g%ranges(4,1) = 0.0_dp
      g%ranges(4,2) = 1.0_dp

      ! SURFACE PRESSURE
      g%ps%name  = 'ps'
      g%ps%id    = NULL_VARID
      g%ps%xtype = NF90_FLOAT
      ! ... dimensions
      g%ps%ndims = 2
      ALLOCATE(g%ps%dim(g%ps%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)
      g%ps%dim(1) = g%lonm%dim(1)
      g%ps%dim(2) = g%latm%dim(1)
      ! ... data
      CALL INIT_NARRAY(g%ps%dat, g%ps%ndims                  &
           ,(/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
           ,VTYPE_REAL)
      IF (ASSOCIATED(aps)) THEN
         DO i=1, g%ps%dim(1)%len
            DO j=1, g%ps%dim(2)%len
               g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/)  &
                    ,(/i,j/)))                           &
                    !qqq
                    !           = aps(i,j)
                    = 101325.0
            END DO
         END DO
      ELSE
         DO i=1, g%ps%dim(1)%len
            DO j=1, g%ps%dim(2)%len
               g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/)  &
                    ,(/i,j/)))                           &
                    = 101325.0
            END DO
         END DO
      END IF
      ! ... attributes
      CALL ADD_NCATT(g%ps, 'long_name', vs='surface pressure')
      CALL ADD_NCATT(g%ps, 'units', vs='Pa')

      ! REFERENCE PRESSURE ...
      g%p0%name  = 'p0'
      g%p0%id    = NULL_VARID
      g%p0%xtype = NF90_FLOAT
      ! ... dimensions
      g%p0%ndims = 0
      ! ... data
      CALL INIT_NARRAY(g%p0%dat, 1, (/ 1 /), VTYPE_REAL)
      g%p0%dat%vr(1) = 1.0
      ! ... attributes
      CALL ADD_NCATT(g%p0, 'long_name', vs='reference pressure')
      CALL ADD_NCATT(g%p0, 'units', vs='Pa')

      ! TIME (MID) ...
      g%timem%name  = 'time'
      g%timem%id    = NULL_VARID
      g%timem%xtype = NF90_FLOAT
      ! ... dimensions
      g%timem%ndims = 1
      ALLOCATE(g%timem%dim(g%timem%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
      g%timem%dim(1)%name  = 'time'
      g%timem%dim(1)%id    = NULL_DIMID
      g%timem%dim(1)%len   = 1
      g%timem%dim(1)%fuid  = .true.
      g%timem%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%timem%dat, g%timem%ndims, (/g%timem%dim(1)%len/) &
           ,VTYPE_DOUBLE)
      ! TIME: SECONDS SINCE MODEL START
      CALL TIME_SPAN_D(dts   &
           , YEAR_START, MONTH_START, DAY_START &
           , HOUR_START, MINUTE_START, SECOND_START  &
           , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
      g%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds
      ! ... attributes
      CALL ADD_NCATT(g%timem, 'long_name'                              &
           ,vs='time in seconds since model start')
      WRITE(tunit, &
           '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
           YEAR_START, MONTH_START, DAY_START &
           , HOUR_START, MINUTE_START, SECOND_START
      CALL ADD_NCATT(g%timem, 'units', vs=TRIM(tunit))

      ! CALCULATE INTs from MIDs
      CALL COMPLETE_GEOHYBGRID(g)

      ! INTERFACE OK
      ok = .true.

    END SUBROUTINE DEFINE_BASEMODEL_GEOHYBGRID
    ! -------------------------------------------------------------------------
#endif

#ifdef COSMO
    ! -------------------------------------------------------------------------
    SUBROUTINE DEFINE_BASEMODEL_GEOHYBGRID(g, ok)

      USE messy_main_mpi_bi,       ONLY: isubpos, my_cart_id, nboundlines
      USE messy_main_data_bi,      ONLY: aps
      USE messy_main_grid_def_bi,  ONLY: hyai, hybi
      USE messy_main_grid_def_mem_bi,ONLY: ie, je, nlev => ke_tot       &
                                       , startlon_tot , startlat_tot    &
                                       , dlon, dlat                     &
                                       , polgam, pollon, pollat, vcoord
      USE messy_main_timer,        ONLY: YEAR_START, MONTH_START, DAY_START &
                                       , HOUR_START, MINUTE_START           &
                                       , SECOND_START                       &
                                       , YEAR, MONTH, DAY                   &
                                       , HOUR, MINUTE, SECOND
      USE messy_main_constants_mem, ONLY: sp
      USE messy_main_timer,         ONLY: TIME_SPAN_D
      USE messy_main_grid_trafo,    ONLY: RGEMPTY, COMPLETE_GEOHYBGRID   &
                                        , PHIROT2PHI, RLAROT2RLA, RGMAX
      USE messy_main_grid_netcdf,   ONLY: ERRMSG, INIT_NCVAR                  &
                                        , NULL_DIMID, NULL_VARID              &
                                        , VTYPE_REAL, VTYPE_DOUBLE &!,VTYPE_INT&
                                        , NF90_FLOAT, NF90_DOUBLE             &
                                        , POSITION, INIT_NARRAY, ADD_NCATT
      USE messy_main_grid,          ONLY: t_geohybgrid, INIT_GEOHYBGRID

      IMPLICIT NONE

      ! I/O
      TYPE(t_geohybgrid), INTENT(OUT) :: g
      LOGICAL           , INTENT(OUT) :: ok

      ! LOCAL
      REAL(DP)                  :: dts
      INTEGER                   :: i,j, n, itot, jtot
      INTEGER                   :: status
      CHARACTER(LEN=100)        :: tunit
      INTEGER                   :: itmp
      REAL(dp)                  :: clon, clat
      REAL(sp), PARAMETER       :: intfac = 100000._sp
      INTEGER                   :: istartlon, istartlat, idlon, idlat, ix
      REAL(sp)                  :: help2
      REAL(dp)                  :: help1

      istartlon = NINT(startlon_tot * intfac)
      istartlat = NINT(startlat_tot * intfac)
      idlon     = NINT(dlon * intfac)
      idlat     = NINT(dlat * intfac)

      ! INIT
      CALL INIT_GEOHYBGRID(g)

      g%name    = 'COSMO'
      g%file    = ''                           ! Filename
      g%t       = 0                            ! time step
      g%corners = 4                            ! number of corners

      g%clonc = .FALSE. ! Axis is not circular
      ! Curvilinear LONGITUDE (MID) ...
      g%clonm%name  = 'lon'
      g%clonm%id    = NULL_VARID
      g%clonm%xtype = NF90_DOUBLE
      ! ... dimensions
      g%clonm%ndims = 2
      ALLOCATE(g%clonm%dim(g%clonm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%clonm%dim(1)%name  = 'lon'
      g%clonm%dim(1)%id    = NULL_DIMID
      g%clonm%dim(1)%len   = ie
      g%clonm%dim(1)%fuid  = .false.
      g%clonm%dim(1)%varid = NULL_VARID
      g%clonm%dim(2)%name  = 'lat'
      g%clonm%dim(2)%id    = NULL_DIMID
      g%clonm%dim(2)%len   = je
      g%clonm%dim(2)%fuid  = .false.
      g%clonm%dim(2)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%clonm%dat, g%clonm%ndims &
           , (/g%clonm%dim(1)%len,g%clonm%dim(2)%len/) &
           ,VTYPE_DOUBLE)

      ! LATITUDE (MID) ...
      g%clatm%name  = 'lat'
      g%clatm%id    = NULL_VARID
      g%clatm%xtype = NF90_DOUBLE
      ! ... dimensions
      g%clatm%ndims = 2
      ALLOCATE(g%clatm%dim(g%clatm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
      g%clatm%dim(1)%name  = 'lon'
      g%clatm%dim(1)%id    = NULL_DIMID
      g%clatm%dim(1)%len   = ie
      g%clatm%dim(1)%fuid  = .false.
      g%clatm%dim(1)%varid = NULL_VARID
      g%clatm%dim(2)%name  = 'lat'
      g%clatm%dim(2)%id    = NULL_DIMID
      g%clatm%dim(2)%len   = je
      g%clatm%dim(2)%fuid  = .false.
      g%clatm%dim(2)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%clatm%dat, g%clatm%ndims &
           , (/g%clatm%dim(1)%len, g%clatm%dim(2)%len/) &
           ,VTYPE_DOUBLE)
      DO i=1, ie
         itot = isubpos(my_cart_id,1) - nboundlines - 1
         itmp = istartlon + (itot+i-1) * idlon
         clon = REAL(itmp,dp) / intfac
         DO j = 1 , je
            jtot  =  isubpos(my_cart_id,2) - nboundlines - 1
            n = (j-1) * g%clonm%dim(1)%len + i
            itmp = istartlat + (jtot+j-1) * idlat
            clat = REAL(itmp,dp) / intfac
            IF  ( (NINT(pollon) == 0 .AND. NINT(polgam)== 180) .OR.  &
                 (NINT(pollon) == 180 .AND. NINT(polgam) == 0)) THEN
               ! geographical coordinates
               g%clonm%dat%vd(n) = clon
               g%clatm%dat%vd(n) = clat
            ELSE
               help1 =  rlarot2rla(clat, clon, pollat, pollon, polgam)
               help2= REAL(help1,sp)
               g%clonm%dat%vd(n) = REAL(help2,dp)
               help1 =  phirot2phi(clat, clon, pollat, polgam)
               help2= REAL(help1,sp)
               g%clatm%dat%vd(n) = REAL(help2,dp)
            END IF
         END DO
      END DO

      CALL ADD_NCATT(g%clonm, 'long_name', vs='mid longitude')
      CALL ADD_NCATT(g%clonm, 'units',     vs='degrees_east')

      ! ... attributes
      CALL ADD_NCATT(g%clatm, 'long_name', vs='mid latitude')
      CALL ADD_NCATT(g%clatm, 'units', vs='degrees_north')

      ! ----------------------------------------------------------
      ! define interfaces:
      ! Curvilinear LONGITUDE (INTERFACES) ...
      g%cloni%name  = 'lon_I'
      g%cloni%id    = NULL_VARID
      g%cloni%xtype = NF90_DOUBLE
      ! ... dimensions
      g%cloni%ndims = 2
      ALLOCATE(g%cloni%dim(g%cloni%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%cloni%dim(1)%name  = 'lon_I'
      g%cloni%dim(1)%id    = NULL_DIMID
      g%cloni%dim(1)%len   = ie+1
      g%cloni%dim(1)%fuid  = .false.
      g%cloni%dim(1)%varid = NULL_VARID
      g%cloni%dim(2)%name  = 'lat_I'
      g%cloni%dim(2)%id    = NULL_DIMID
      g%cloni%dim(2)%len   = je+1
      g%cloni%dim(2)%fuid  = .false.
      g%cloni%dim(2)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%cloni%dat, g%cloni%ndims &
           , (/g%cloni%dim(1)%len,g%cloni%dim(2)%len/) &
           ,VTYPE_DOUBLE)

      ! LATITUDE (MID) ...
      g%clati%name  = 'lat_I'
      g%clati%id    = NULL_VARID
      g%clati%xtype = NF90_DOUBLE

      ! ... dimensions
      g%clati%ndims = 2
      ALLOCATE(g%clati%dim(g%clati%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
      g%clati%dim(1)%name  = 'lon_I'
      g%clati%dim(1)%id    = NULL_DIMID
      g%clati%dim(1)%len   = ie+1
      g%clati%dim(1)%fuid  = .false.
      g%clati%dim(1)%varid = NULL_VARID
      g%clati%dim(2)%name  = 'lat_I'
      g%clati%dim(2)%id    = NULL_DIMID
      g%clati%dim(2)%len   = je+1
      g%clati%dim(2)%fuid  = .false.
      g%clati%dim(2)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%clati%dat, g%clati%ndims &
           , (/g%clati%dim(1)%len, g%clati%dim(2)%len/) &
           ,VTYPE_DOUBLE)

      DO i=1, ie+1
         itot = isubpos(my_cart_id,1) - nboundlines - 1
         itmp = istartlon + (itot+i) * idlon - NINT(1.5_dp * dlon * intfac)
         clon = REAL(itmp,dp) / intfac
         DO j = 1 , je+1
            jtot  =  isubpos(my_cart_id,2) - nboundlines - 1
            n = (j-1) * g%cloni%dim(1)%len + i
            itmp = istartlat + (jtot+j) * idlat - NINT(1.5_dp * dlat * intfac)
            clat = REAL(itmp,dp) / intfac
            IF  ( (NINT(pollon) == 0 .AND. NINT(polgam)== 180) .OR.  &
                 (NINT(pollon) == 180 .AND. NINT(polgam) == 0)) THEN
               ! geographical coordinates
               g%cloni%dat%vd(n) = clon
               g%clati%dat%vd(n) = clat
            ELSE
               help1 =  rlarot2rla(clat, clon, pollat, pollon, polgam)
               help2= REAL(help1,sp)
               g%cloni%dat%vd(n) = REAL(help2,dp)
               help1 =  phirot2phi(clat, clon, pollat, polgam)
               help2= REAL(help1,sp)
               g%clati%dat%vd(n) = REAL(help2,dp)
            ENDIF
         END DO
      END DO

      CALL ADD_NCATT(g%cloni, 'long_name', vs='interface longitude')
      CALL ADD_NCATT(g%cloni, 'units', vs='degrees_east')

      ! ... attributes
      CALL ADD_NCATT(g%clati, 'long_name', vs='interface latitude')
      CALL ADD_NCATT(g%clati, 'units', vs='degrees_north')

      g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
      g%ranges(1,2) = RGEMPTY

      g%ranges(2,1) = RGEMPTY !-90.0_sp
      g%ranges(2,2) = RGEMPTY

      ! DETERMIN MIN / MAX LONGITUDE
      ! CHECK IF PE overlaps 360:0 degree line, but has gap between [0,360]
      IF (.NOT. ANY (ABS(g%cloni%dat%vd(:)) < 20._dp) .AND. &
           ANY(ABS(g%cloni%dat%vd(:)) > 160._dp) ) THEN
         DO ix = 1, SIZE(g%cloni%dat%vd)
            IF (g%cloni%dat%vd(ix) < 0._dp) THEN
               g%cloni%dat%vd(ix) = g%cloni%dat%vd(ix) +360._dp
            END IF
         END DO
      END IF
      g%minmaxlonlat(1,1) = MINVAL(g%cloni%dat%vd)
      g%minmaxlonlat(1,2) = MAXVAL(g%cloni%dat%vd)

      ! DETERMIN MIN / MAX  LATITUDE
      g%minmaxlonlat(2,1) = MINVAL(g%clati%dat%vd)
      g%minmaxlonlat(2,2) = MAXVAL(g%clati%dat%vd)
      !**************************************************************************
      ! Rotated Pole
      ! ... longitude
      CALL INIT_NCVAR(g%pollon)
      g%pollon%name  = TRIM('pollon')
      g%pollon%xtype = NF90_FLOAT
      !
      ! NO DIMENSIONS !!!
      g%pollon%ndims = 1
      ALLOCATE(g%pollon%dim(1), STAT=status)
      CALL ERRMSG(substr,status,4)
      g%pollon%dim(1)%name  = TRIM(g%pollon%name)//'_dim'
      g%pollon%dim(1)%id    = NULL_DIMID
      g%pollon%dim(1)%len   = 1
      g%pollon%dim(1)%fuid  = .false.
      g%pollon%dim(1)%varid = NULL_VARID
      !
      ! ATTRIBUTES ...
      ! ... LONGNAME ATTRIBUTE
      CALL ADD_NCATT(g%pollon, 'long_name', vs='g_north_pole_longitude')
      ! ... UNITS ATTRIBUT
      CALL ADD_NCATT(g%pollon, 'units' ,vs='degree_east')
      !
      ! DATA
      CALL INIT_NARRAY(g%pollon%dat, 1, (/ 1 /), VTYPE_REAL)
      g%pollon%dat%vr(1) = REAL(pollon,sp)
      ! ... latitude
      CALL INIT_NCVAR(g%pollat)
      g%pollat%name  = TRIM('pollat')
      g%pollat%xtype = NF90_FLOAT
      !
      ! NO DIMENSIONS !!!
      g%pollat%ndims = 1
      ALLOCATE(g%pollat%dim(1), STAT=status)
      CALL ERRMSG(substr,status,4)
      g%pollat%dim(1)%name  = TRIM(g%pollat%name)//'_dim'
      g%pollat%dim(1)%id    = NULL_DIMID
      g%pollat%dim(1)%len   = 1
      g%pollat%dim(1)%fuid  = .false.
      g%pollat%dim(1)%varid = NULL_VARID
      !
      ! ATTRIBUTES ...
      ! ... LONGNAME ATTRIBUTE
      CALL ADD_NCATT(g%pollat, 'long_name', vs='g_north_pole_latitude')
      ! ... UNITS ATTRIBUT
      CALL ADD_NCATT(g%pollat, 'units' ,vs='degree_north')
      !
      ! DATA
      CALL INIT_NARRAY(g%pollat%dat, 1, (/ 1 /), VTYPE_REAL)
      g%pollat%dat%vr(1) = REAL(pollat,sp)
      ! ... angle
      CALL INIT_NCVAR(g%polgam)
      g%polgam%name  = TRIM('polgam')
      g%polgam%xtype = NF90_FLOAT
      !
      ! NO DIMENSIONS !!!
      g%polgam%ndims = 1
      ALLOCATE(g%polgam%dim(1), STAT=status)
      CALL ERRMSG(substr,status,4)
      g%polgam%dim(1)%name  = TRIM(g%polgam%name)//'_dim'
      g%polgam%dim(1)%id    = NULL_DIMID
      g%polgam%dim(1)%len   = 1
      g%polgam%dim(1)%fuid  = .false.
      g%polgam%dim(1)%varid = NULL_VARID
      !
      ! ATTRIBUTES ...
      ! ... LONGNAME ATTRIBUTE
      CALL ADD_NCATT(g%polgam, 'long_name', vs='g_north_pole_angle')
      ! ... UNITS ATTRIBUT
      CALL ADD_NCATT(g%polgam, 'units' ,vs='-')
      !
      ! DATA
      CALL INIT_NARRAY(g%polgam%dat, 1, (/ 1 /), VTYPE_REAL)
      g%polgam%dat%vr(1) = REAL(polgam,sp)

      ! Rotated  Curvilinear LONGITUDE (MID) ...
      g%rlonm%name  = 'lon'
      g%rlonm%id    = NULL_VARID
      g%rlonm%xtype = NF90_FLOAT
      g%rlonc = .FALSE. ! Axis is not circular
      ! ... dimensions
      g%rlonm%ndims = 1
      ALLOCATE(g%rlonm%dim(g%rlonm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%rlonm%dim(1)%name  = 'lon'
      g%rlonm%dim(1)%id    = NULL_DIMID
      g%rlonm%dim(1)%len   = ie
      g%rlonm%dim(1)%fuid  = .false.

      ! ... data
      CALL INIT_NARRAY(g%rlonm%dat, g%rlonm%ndims &
           , (/g%rlonm%dim(1)%len/),VTYPE_REAL)
      DO i=1, ie
         itot = isubpos(my_cart_id,1) - nboundlines - 1
         itmp = istartlon + (itot+i-1) * idlon
         g%rlonm%dat%vr(i) = REAL(itmp,sp) / intfac
      END DO

      ! ... attributes
      CALL ADD_NCATT(g%rlonm, 'long_name', vs='longitude')
      CALL ADD_NCATT(g%rlonm, 'units', vs='degrees_east')
      g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
      g%ranges(1,2) = RGEMPTY

      ! LATITUDE (MID) ...
      g%rlatm%name  = 'lat'
      g%rlatm%id    = NULL_VARID
      g%rlatm%xtype = NF90_FLOAT
      ! ... dimensions
      g%rlatm%ndims = 1
      ALLOCATE(g%rlatm%dim(g%rlatm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
      g%rlatm%dim(1)%name  = 'lat'
      g%rlatm%dim(1)%id    = NULL_DIMID
      g%rlatm%dim(1)%len   = je
      g%rlatm%dim(1)%fuid  = .false.
      g%rlatm%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%rlatm%dat, g%rlatm%ndims &
           , (/g%rlatm%dim(1)%len/),VTYPE_REAL)
      DO j=1,je
         jtot = isubpos(my_cart_id,2) - nboundlines - 1
         itmp = istartlat + (jtot+j-1) * idlat
         g%rlatm%dat%vr(j) = REAL(itmp,sp) / intfac
      END DO

      ! ... attributes
      CALL ADD_NCATT(g%rlatm, 'long_name', vs='mid latitude')
      CALL ADD_NCATT(g%rlatm, 'units', vs='degrees_north')

      ! ----------------------------------------------------------
      ! define interfaces:
      ! Curvilinear LONGITUDE (INTERFACES) ...
      g%rloni%name  = 'lon_I'
      g%rloni%id    = NULL_VARID
      g%rloni%xtype = NF90_FLOAT
      ! ... dimensions
      g%rloni%ndims = 1
      ALLOCATE(g%rloni%dim(g%rloni%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%rloni%dim(1)%name  = 'lon_I'
      g%rloni%dim(1)%id    = NULL_DIMID
      g%rloni%dim(1)%len   = ie+1
      g%rloni%dim(1)%fuid  = .false.
      g%rloni%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%rloni%dat, g%rloni%ndims &
           , (/g%rloni%dim(1)%len/)  ,VTYPE_REAL)

      DO i=1, ie+1
         itot = isubpos(my_cart_id,1) - nboundlines - 1
         itmp = istartlon + (itot+i) * idlon - NINT(1.5_dp * dlon * intfac)
         g%rloni%dat%vr(i) = REAL(itmp,sp) / intfac
      END DO

      ! ... attributes
      CALL ADD_NCATT(g%rloni, 'long_name', vs='interface longitude')
      CALL ADD_NCATT(g%rloni, 'units', vs='degrees_east')

      ! LATITUDE (MID) ...
      g%rlati%name  = 'lat_I'
      g%rlati%id    = NULL_VARID
      g%rlati%xtype = NF90_FLOAT
      ! ... dimensions
      g%rlati%ndims = 1
      ALLOCATE(g%rlati%dim(g%rlati%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
      g%rlati%dim(1)%name  = 'lat_I'
      g%rlati%dim(1)%id    = NULL_DIMID
      g%rlati%dim(1)%len   = je+1
      g%rlati%dim(1)%fuid  = .false.
      g%rlati%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%rlati%dat, g%rlati%ndims &
           , (/g%rlati%dim(1)%len/) ,VTYPE_REAL)
      DO j = 1 , je+1
         jtot  =  isubpos(my_cart_id,2) - nboundlines - 1
         n = (j-1) * g%cloni%dim(1)%len + i
         itmp = istartlat + (jtot+j) * idlat - NINT(1.5_dp * dlat * intfac)
         g%rlati%dat%vr(j) = REAL(itmp,sp) / intfac
      END DO

      ! ... attributes
      CALL ADD_NCATT(g%rlati, 'long_name', vs='interface latitude')
      CALL ADD_NCATT(g%rlati, 'units', vs='degrees_north')

      ! **********************************************************************
      IF (vcoord%ivctype == 1) THEN
         write (*,*) 'DEFINING BASEGRID as PRESSURE GRID'

         ! HYBRID-A-COEFFICIENTS (INTERFACES) ...
         g%hyai%name  = 'hyai'
         g%hyai%id    = NULL_VARID
         g%hyai%xtype = NF90_DOUBLE
         ! ... dimensions
         g%hyai%ndims = 1
         ALLOCATE(g%hyai%dim(g%hyai%ndims), STAT=status)
         CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
         g%hyai%dim(1)%name  = 'ilev'
         g%hyai%dim(1)%id    = NULL_DIMID
         g%hyai%dim(1)%len   = nlev+1
         g%hyai%dim(1)%fuid  = .false.
         g%hyai%dim(1)%varid = NULL_VARID
         ! ... data
         CALL INIT_NARRAY(g%hyai%dat, g%hyai%ndims, (/g%hyai%dim(1)%len/) &
              ,VTYPE_DOUBLE)
         g%hyai%dat%vd(:) = hyai
         ! ... attributes
         CALL ADD_NCATT(g%hyai, 'long_name'                              &
              ,vs='hybrid-A-coefficients at layer interfaces')
         CALL ADD_NCATT(g%hyai, 'units', vs='1')
         g%ranges(3,1) = 0.0_sp
         g%ranges(3,2) = 0.0_sp

         ! ********************************************************************
         ! HYBRID-B-COEFFICIENTS (INTERFACES) ...
         g%hybi%name  = 'hybi'
         g%hybi%id    = NULL_VARID
         g%hybi%xtype = NF90_DOUBLE
         ! ... dimensions
         g%hybi%ndims = 1
         ALLOCATE(g%hybi%dim(g%hybi%ndims), STAT=status)
         CALL ERRMSG('INTERFACE_GEOHYBGRID',status,4)
         g%hybi%dim(1) = g%hyai%dim(1)
         ! ... data
         CALL INIT_NARRAY(g%hybi%dat, g%hybi%ndims, (/g%hybi%dim(1)%len/) &
              ,VTYPE_DOUBLE)
         g%hybi%dat%vd(:) = hybi !vct(nvclev+1:2*nvclev)
         ! ... attributes
         CALL ADD_NCATT(g%hybi, 'long_name'                              &
              ,vs='hybrid-B-coefficients at layer interfaces')
         CALL ADD_NCATT(g%hybi, 'units', vs='1')
         g%ranges(4,1) = 0.0_dp
         g%ranges(4,2) = 1.0_dp

      ELSE

         !write (*,*) 'DEFINING BASEGRID as HEIGHT GRID'

         ! prepare space for the up-to-the-minte 3D pressure field
         ! required for regridding of pressure coordinates to
         ! height coordinates
         ! 3d PRESSURE FIELD INTERFACES
         g%pressi%name  = 'pressi'
         g%pressi%id    = NULL_VARID
         g%pressi%xtype = NF90_DOUBLE
         ! ... dimensions
         g%pressi%ndims = 3
         ALLOCATE(g%pressi%dim(g%pressi%ndims), STAT=status)
         CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)

         g%pressi%dim(1) = g%clonm%dim(1)
         g%pressi%dim(2) = g%clonm%dim(2)
         g%pressi%dim(3)%name  = 'ilev'
         g%pressi%dim(3)%id    = NULL_DIMID
         g%pressi%dim(3)%len   = nlev+1
         g%pressi%dim(3)%fuid  = .false.
         g%pressi%dim(3)%varid = NULL_VARID

         ! ... data
         CALL INIT_NARRAY(g%pressi%dat, g%pressi%ndims      &
              ,(/g%pressi%dim(1)%len, g%pressi%dim(2)%len   &
              ,  g%pressi%dim(3)%len/)                      &
              , VTYPE_DOUBLE)
         g%ranges(3,1) = 0._dp
         ! define the possible maximum (1085 hPa in the Siberian low)
         g%ranges(3,2) = RGMAX

         ! 3d PRESSURE FIELD MIDS
         g%pressm%name  = 'pressm'
         g%pressm%id    = NULL_VARID
         g%pressm%xtype = NF90_DOUBLE
         ! ... dimensions
         g%pressm%ndims = 3
         ALLOCATE(g%pressm%dim(g%pressm%ndims), STAT=status)
         CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)

         g%pressm%dim(1) = g%clonm%dim(1)
         g%pressm%dim(2) = g%clonm%dim(2)
         g%pressm%dim(3)%name  = 'ilev'
         g%pressm%dim(3)%id    = NULL_DIMID
         g%pressm%dim(3)%len   = nlev
         g%pressm%dim(3)%fuid  = .false.
         g%pressm%dim(3)%varid = NULL_VARID

         ! ... data
         CALL INIT_NARRAY(g%pressm%dat, g%pressm%ndims      &
              ,(/g%pressm%dim(1)%len, g%pressm%dim(2)%len   &
              ,  g%pressm%dim(3)%len/)                      &
              , VTYPE_DOUBLE)

      END IF

      ! ***********************************************************************
      ! SURFACE PRESSURE
      g%ps%name  = 'ps'
      g%ps%id    = NULL_VARID
      g%ps%xtype = NF90_FLOAT
      ! ... dimensions
      g%ps%ndims = 2
      ALLOCATE(g%ps%dim(g%ps%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)

      g%ps%dim(1) = g%clonm%dim(1)
      g%ps%dim(2) = g%clonm%dim(2)

      ! ... data
      CALL INIT_NARRAY(g%ps%dat, g%ps%ndims                  &
           ,(/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
           ,VTYPE_REAL)
      IF (ASSOCIATED(aps)) THEN
         DO i=1, g%ps%dim(1)%len
            DO j=1, g%ps%dim(2)%len
               g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/)  &
                    ,(/i,j/)))                           &
                    !qqq
                    ! = aps(i,j)
                    = 101325.0_sp
            END DO
         END DO
      ELSE
         DO i=1, g%ps%dim(1)%len
            DO j=1, g%ps%dim(2)%len
               g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/)  &
                    ,(/i,j/)))                           &
                    = 101325.0_sp
            END DO
         END DO
      END IF
      ! ... attributes
      CALL ADD_NCATT(g%ps, 'long_name', vs='surface pressure')
      CALL ADD_NCATT(g%ps, 'units', vs='Pa')

      ! ***********************************************************************
      ! REFERENCE PRESSURE ...
      g%p0%name  = 'p0'
      g%p0%id    = NULL_VARID
      g%p0%xtype = NF90_FLOAT
      ! ... dimensions
      g%p0%ndims = 0
      ! ... data
      CALL INIT_NARRAY(g%p0%dat, 1, (/ 1 /), VTYPE_REAL)
      g%p0%dat%vr(1) = 1.0_sp
      ! ... attributes
      CALL ADD_NCATT(g%p0, 'long_name', vs='reference pressure')
      CALL ADD_NCATT(g%p0, 'units', vs='Pa')

      ! ***********************************************************************
!!$      ! DEFINE MASK of defined points
!!$      g%imask%name = 'mask'
!!$      g%imask%id    = NULL_VARID
!!$      g%imask%xtype = NF90_DOUBLE
!!$      ! ... dimensions
!!$      g%imask%ndims = 2
!!$      ALLOCATE(g%imask%dim(g%imask%ndims), STAT=status)
!!$      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
!!$      g%imask%dim(1)%name = 'lon'
!!$      g%imask%dim(1)%id    = NULL_DIMID
!!$      g%imask%dim(1)%len   = ie
!!$      g%imask%dim(1)%fuid  = .false.
!!$      g%imask%dim(1)%varid = NULL_VARID
!!$      g%imask%dim(2)%name  = 'lat'
!!$      g%imask%dim(2)%id    = NULL_DIMID
!!$      g%imask%dim(2)%len   = je
!!$      g%imask%dim(2)%fuid  = .false.
!!$      g%imask%dim(2)%varid = NULL_VARID
!!$
!!$      ! .. ATTRIBUTES
!!$      CALL ADD_NCATT(g%imask, 'long_name', vs='mask of defined points')
!!$
!!$      ! ... data
!!$      CALL INIT_NARRAY(g%imask%dat, g%imask%ndims &
!!$           , (/g%imask%dim(1)%len,g%imask%dim(2)%len/) &
!!$           ,VTYPE_INT)
!!$
!!$      g%imask%dat%vi = 1

      ! ***********************************************************************
      ! TIME (MID) ...
      g%timem%name  = 'time'
      g%timem%id    = NULL_VARID
      g%timem%xtype = NF90_DOUBLE
      ! ... dimensions
      g%timem%ndims = 1
      ALLOCATE(g%timem%dim(g%timem%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
      g%timem%dim(1)%name  = 'time'
      g%timem%dim(1)%id    = NULL_DIMID
      g%timem%dim(1)%len   = 1
      g%timem%dim(1)%fuid  = .true.
      g%timem%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%timem%dat, g%timem%ndims, (/g%timem%dim(1)%len/) &
           ,VTYPE_DOUBLE)
      ! TIME: SECONDS SINCE MODEL START
      CALL TIME_SPAN_D(dts   &
           , YEAR_START, MONTH_START, DAY_START &
           , HOUR_START, MINUTE_START, SECOND_START  &
           , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
      g%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds

      ! ... attributes
      CALL ADD_NCATT(g%timem, 'long_name'                              &
           ,vs='time in seconds since model start')
      WRITE(tunit, &
           '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
           YEAR_START, MONTH_START, DAY_START &
           , HOUR_START, MINUTE_START, SECOND_START

      CALL ADD_NCATT(g%timem, 'units', vs=TRIM(tunit))

      ! CALCULATE INTs from MIDs
      CALL COMPLETE_GEOHYBGRID(g)

      ! INTERFACE OK
      ok = .true.

    END SUBROUTINE DEFINE_BASEMODEL_GEOHYBGRID
    ! -------------------------------------------------------------------------
#endif
!COSMO

#ifdef CESM1
    ! -------------------------------------------------------------------------
    SUBROUTINE DEFINE_BASEMODEL_GEOHYBGRID(g, ok)

      USE messy_main_mpi_bi,          ONLY: dcg, SCATTER_GP
      USE messy_main_grid_def_mem_bi, ONLY: nvclev, vct, nlev &
#ifndef HOMMESE
                                          , nlon, nlat &
#endif
                                          , nproma, ngpblks, pcols, npromz
      USE messy_main_data_bi,       ONLY: aps
      USE filenames,                ONLY: ncdata
      USE messy_main_grid_def_bi,   ONLY: philat_2d, philon_2d
      USE messy_main_timer,         ONLY: YEAR_START, MONTH_START, DAY_START &
                                        , HOUR_START, MINUTE_START &
                                        , SECOND_START             &
                                        , YEAR, MONTH, DAY         &
                                        , HOUR, MINUTE, SECOND
      USE messy_main_timer,         ONLY: TIME_SPAN_D
      USE messy_main_constants_mem, ONLY: api=>pi, sp
      USE messy_main_grid_trafo,    ONLY: RGEMPTY, COMPLETE_GEOHYBGRID
      USE messy_main_grid_tools,    ONLY: DEALINE_ARRAY &
                                        , RGTOOL_CONVERT_DAT2PREDEFVAR
      USE messy_main_grid_netcdf,   ONLY: ERRMSG  &
                                        , NULL_DIMID, NULL_VARID &
                                        , VTYPE_REAL, VTYPE_DOUBLE, VTYPE_INT &
                                        , nf90_float, NF90_DOUBLE, NF90_INT  &
                                        , POSITION, INIT_NARRAY, ADD_NCATT &
                                        , INIT_NCVAR, IMPORT_NCVAR         &
                                        , t_ncvar, PRINT_NCVAR, COPY_NCATT &
                                        , double_narray
      USE messy_main_grid,          ONLY: t_geohybgrid, INIT_GEOHYBGRID &
                                        , PRINT_GEOHYBGRID

      IMPLICIT NONE

      ! I/O
      TYPE (t_geohybgrid), INTENT(OUT) :: g
      LOGICAL            , INTENT(OUT) :: ok

      ! LOCAL
      REAL(DP) :: dts
      INTEGER :: i,j, n, ix,jx,  dim1, gsize
      INTEGER :: status
      CHARACTER(LEN=100)       :: tunit
      ! IMPORTED LONGITUDES / LATITUDES OF CENTERS / CORNERS
      TYPE(t_ncvar)            :: loclon, loclat
      ! LOCAL 2D LONGITUDE / LATITUDE FIELDS
      REAL(DP), DIMENSION(:,:), POINTER :: loclon_2d => NULL()
      REAL(DP), DIMENSION(:,:), POINTER :: loclat_2d => NULL()

      ! POINTER FOR FIELD TRANSFORMATION
      REAL(DP), DIMENSION(:,:),     POINTER :: dummy2d  => NULL()
      REAL(DP), DIMENSION(:,:,:,:), POINTER :: dum2dloc => NULL()

      ! POINTER FOR SCATTER ROUTINES
      REAL(DP), DIMENSION(:),   POINTER   :: gdat   => NULL()
      REAL(DP), DIMENSION(:,:), POINTER   :: ldat   => NULL()
      REAL(dp), DIMENSION(:), ALLOCATABLE :: test
      REAL(dp)                            :: diff1, diff2
      REAL(dp)                            :: minl, maxl

      ! INIT
      CALL INIT_GEOHYBGRID(g)

      g%name    = 'CESM1'
      g%file    = ''                              ! Filename
      g%t       = 0                               ! time step
      g%corners = 5                               ! number of corners

#ifndef HOMMESE
      ! LONGITUDE (MID) ...
      g%lonm%name  = 'lon'
      g%lonm%id    = NULL_VARID
      g%lonm%xtype = NF90_FLOAT
      ! ... dimensions
      g%lonm%ndims = 1
      ALLOCATE(g%lonm%dim(g%lonm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%lonm%dim(1)%name  = 'lon'
      g%lonm%dim(1)%id    = NULL_DIMID
      g%lonm%dim(1)%len   = nlon
      g%lonm%dim(1)%fuid  = .false.
      g%lonm%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%lonm%dat, g%lonm%ndims, (/g%lonm%dim(1)%len/) &
           ,VTYPE_REAL)
      DO i=1, nlon
         g%lonm%dat%vr(i) = 360. * REAL(i-1) / REAL(nlon)
      END DO
      ! ... attributes
      CALL ADD_NCATT(g%lonm, 'long_name', vs='longitude')
      CALL ADD_NCATT(g%lonm, 'units', vs='degrees_east')
      g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
      g%ranges(1,2) = RGEMPTY

      ! LATITUDE (MID) ...
      g%latm%name  = 'lat'
      g%latm%id    = NULL_VARID
      g%latm%xtype = NF90_FLOAT
      ! ... dimensions
      g%latm%ndims = 1
      ALLOCATE(g%latm%dim(g%latm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
      g%latm%dim(1)%name  = 'lat'
      g%latm%dim(1)%id    = NULL_DIMID
      g%latm%dim(1)%len   = nlat
      g%latm%dim(1)%fuid  = .false.
      g%latm%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%latm%dat, g%latm%ndims, (/g%latm%dim(1)%len/) &
           ,VTYPE_REAL)
      DO j=1,nlat
         g%latm%dat%vr(j)=philat(j)!-ABS(philat(1)-philat(2))
         !     g%latm%dat%vr(ngl+1-j)=-g%latm%dat%vr(j)
      END DO
      ! ... attributes
      CALL ADD_NCATT(g%latm, 'long_name', vs='latitude')
      CALL ADD_NCATT(g%latm, 'units', vs='degrees_north')
      g%ranges(2,1) = -90.0_dp
      g%ranges(2,2) =  90.0_dp

#else

      ! unstructured grid
      g%lonc = .FALSE. ! Axis is not circular

      ! DO not define IMASK for CESM: leads to decomposition dependent
      ! errors in IMPORT_GRID

      ! Unstructured LONGITUDE Centers(MID) ...
      CALL INIT_NCVAR(g%ulonm)
      CALL INIT_NCVAR(g%ulatm)

      g%ulonm%name  = 'unstruc_lon'
      g%ulonm%id    = NULL_VARID
      g%ulonm%xtype = NF90_DOUBLE
      ! ... dimensions
      g%ulonm%ndims = 2
      ALLOCATE(g%ulonm%dim(g%ulonm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%ulonm%dim(1)%name = 'unstruc_lon'
      g%ulonm%dim(1)%id    = NULL_DIMID
      g%ulonm%dim(1)%len   = SIZE(philon_2d,1)
      g%ulonm%dim(1)%fuid  = .false.
      g%ulonm%dim(1)%varid = NULL_VARID
      g%ulonm%dim(2)%name  = 'unstruc_lat'
      g%ulonm%dim(2)%id    = NULL_DIMID
      g%ulonm%dim(2)%len   = SIZE(philon_2d,2)
      g%ulonm%dim(2)%fuid  = .false.
      g%ulonm%dim(2)%varid = NULL_VARID

      g%ranges(1,1) = RGEMPTY
      g%ranges(1,2) = RGEMPTY
      ! .. ATTRIBUTES
      CALL ADD_NCATT(g%ulonm, 'long_name', vs='longitude')
      CALL ADD_NCATT(g%ulonm, 'units', vs='degrees_east')

      ! ... data
      CALL INIT_NARRAY(g%ulonm%dat, g%ulonm%ndims &
           , (/g%ulonm%dim(1)%len,g%ulonm%dim(2)%len/) &
           ,VTYPE_DOUBLE)

      CALL INIT_NCVAR(loclon)
      CALL IMPORT_NCVAR(loclon,varname='grid_center_lon', file=ncdata)

      CALL DOUBLE_NARRAY(loclon%dat)

      ALLOCATE(loclon_2d(nproma, ngpblks))
      ! This scatters from global 1D  to local 2D !
      CALL SCATTER_GP(loclon%dat%vd,loclon_2d,dcg)

      ! FLAG empty points
      DO jx = 1, ngpblks
         IF (npromz(jx) /= nproma) THEN
            DO ix = npromz(jx)+1, nproma
               loclon_2d(ix,jx) = loclon_2d(npromz(jx),jx)
            END DO
         END IF
      END DO
      ! ... assign 1d GRID Data
      dim1=g%ulonm%dim(1)%len
      DO ix = 1, g%ulonm%dim(1)%len
         DO jx = 1, g%ulonm%dim(2)%len
            g%ulonm%dat%vd((jx-1)*dim1 + ix) = loclon_2d(ix,jx)
         END DO
      END DO

      CALL INIT_NCVAR(loclon)
      DEALLOCATE(loclon_2d)
      NULLIFY(loclon_2d)

      ! -------------------------------------
      ! Unstructured LATITUDE Centers(MID) ...
      g%ulatm%name  = 'unstruc_lat'
      g%ulatm%id    = NULL_VARID
      g%ulatm%xtype = NF90_DOUBLE
      ! ... dimensions
      g%ulatm%ndims = 2
      ALLOCATE(g%ulatm%dim(g%ulatm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%ulatm%dim(1)%name = 'unstruc_lon'
      g%ulatm%dim(1)%id    = NULL_DIMID
      g%ulatm%dim(1)%len   = SIZE(philat_2d,1)
      g%ulatm%dim(1)%fuid  = .false.
      g%ulatm%dim(1)%varid = NULL_VARID
      g%ulatm%dim(2)%name  = 'unstruc_lat'
      g%ulatm%dim(2)%id    = NULL_DIMID
      g%ulatm%dim(2)%len   = SIZE(philat_2d,2)
      g%ulatm%dim(2)%fuid  = .false.
      g%ulatm%dim(2)%varid = NULL_VARID

      g%ranges(2,1) = RGEMPTY !-90.0_dp
      g%ranges(2,2) = RGEMPTY ! 90.0_dp
      ! .. ATTRIBUTES
      CALL ADD_NCATT(g%ulatm, 'long_name', vs='latitude')
      CALL ADD_NCATT(g%ulatm, 'units', vs='degrees_north')

      ! ... data
      CALL INIT_NARRAY(g%ulatm%dat, g%ulatm%ndims &
           , (/g%ulatm%dim(1)%len,g%ulatm%dim(2)%len/) &
           ,VTYPE_DOUBLE)

      CALL INIT_NCVAR(loclat)
      CALL IMPORT_NCVAR(loclat,varname='grid_center_lat', file=ncdata)
      CALL DOUBLE_NARRAY(loclat%dat)
      ALLOCATE(loclat_2d(nproma, ngpblks))
      ! This scatters from global 1D  to local 2D !
      CALL SCATTER_GP(loclat%dat%vd,loclat_2d,dcg)

      ! FLAG empty points
      DO jx = 1, ngpblks
         IF (npromz(jx) /= nproma) THEN
            DO ix = npromz(jx)+1, nproma
               loclat_2d(ix,jx) = loclat_2d(npromz(jx),jx)
            END DO
         END IF
      END DO
      ! ... assign 1d GRID Dat
      dim1=g%ulatm%dim(1)%len
      DO ix = 1, g%ulatm%dim(1)%len
         DO jx = 1, g%ulatm%dim(2)%len
            g%ulatm%dat%vd((jx-1)*dim1 + ix) = loclat_2d(ix,jx)
         END DO
      END DO

      CALL INIT_NCVAR(loclat)
      DEALLOCATE(loclat_2d)
      NULLIFY(loclat_2d)
      ! -------------------------------------

      ! -------------------------------------
      CALL INIT_NCVAR(g%uloni)
      CALL INIT_NCVAR(g%ulati)

      ! -------------------------------------
      ! -------------------------------------
      ! Unstructured LONGITUDE corners (INTERFACES) ...
      g%uloni%name  = 'unstruc_corners_lon'
      g%uloni%id    = NULL_VARID
      g%uloni%xtype = NF90_DOUBLE
      ! ... dimensions
      g%uloni%ndims = 3
      ALLOCATE(g%uloni%dim(g%uloni%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%uloni%dim(1)%name = 'unstruc_lon'
      g%uloni%dim(1)%id    = NULL_DIMID
      g%uloni%dim(1)%len   = g%ulonm%dim(1)%len
      g%uloni%dim(1)%fuid  = .false.
      g%uloni%dim(1)%varid = NULL_VARID
      g%uloni%dim(2)%name = 'unstruc_lon'
      g%uloni%dim(2)%id    = NULL_DIMID
      g%uloni%dim(2)%len   = g%ulonm%dim(2)%len
      g%uloni%dim(2)%fuid  = .false.
      g%uloni%dim(2)%varid = NULL_VARID
      g%uloni%dim(3)%name  = 'num_corners'
      g%uloni%dim(3)%id    = NULL_DIMID
      g%uloni%dim(3)%len   = g%corners ! or loclon%dim(1)%len
      g%uloni%dim(3)%fuid  = .false.
      g%uloni%dim(3)%varid = NULL_VARID

      ! COPY ATTRIBUTES
      CALL ADD_NCATT(g%uloni, 'long_name', vs='corner longitude')
      CALL ADD_NCATT(g%uloni, 'units',     vs='degrees_east')

      ! ... data
      CALL INIT_NARRAY(g%uloni%dat, g%uloni%ndims &
           , (/g%uloni%dim(1)%len,g%uloni%dim(2)%len,g%uloni%dim(3)%len/) &
           ,VTYPE_DOUBLE)

      CALL INIT_NCVAR(loclon)
      CALL IMPORT_NCVAR(loclon,varname='grid_corner_lon', file=ncdata)
      CALL DOUBLE_NARRAY(loclon%dat)
      ! convert from global 1D to  global (corners,ncol)
      CALL DEALINE_ARRAY(status,loclon%dim(1)%len, loclon%dim(2)%len &
           , loclon%dat%vd, dummy2d)

      ! convert from  global (corners,ncol) to local (nproma, ngpblks, corners)
      ALLOCATE(dum2dloc(nproma,ngpblks,5,1))
      DO ix = 1,SIZE(dummy2d,1)
         gdat => dummy2d(ix,:)
         ldat => dum2dloc(:,:,ix,1)
         CALL SCATTER_GP(gdat,ldat,dcg)
      END DO

      ! FLAG empty points
      DO jx = 1, ngpblks
         IF (npromz(jx) /= nproma) THEN
            DO ix = npromz(jx)+1, nproma
               dum2dloc(ix,jx,:,:) = dum2dloc(npromz(jx),jx,:,:)
            END DO
         END IF
      END DO

      ! ... assign 1d GRID Data
      dim1=g%uloni%dim(1)%len
      gsize=g%uloni%dim(1)%len *  g%uloni%dim(2)%len
      DO n = 1, g%corners
         DO ix = 1, g%uloni%dim(1)%len
            DO jx = 1, g%uloni%dim(2)%len
               g%uloni%dat%vd((jx-1)*dim1 + ix + (n-1)*gsize) = &
                    dum2dloc(ix,jx,n,1)
            END DO
         END DO
      END DO

      ! SET MINIMUM and MAXIMUM for GRID reduction of input grids
      g%minmaxlonlat(1,1) = MINVAL(dum2dloc)
      g%minmaxlonlat(1,2) = MAXVAL(dum2dloc)

      ! CHECK IF PE overlaps 360:0 degree line, but has gap between [0,360]
      ! COULD BE THE case IF DIFF > 270 ....)
      diff1 = g%minmaxlonlat(1,2)-g%minmaxlonlat(1,1)

      IF (diff1 > 180._dp) then
         ! CHECK, if PE covers domain crossing the 0 degree line
         ! TODO: is there a better way?
         ALLOCATE(test(size(g%uloni%dat%vd)))
         test(:) = g%uloni%dat%vd(:)
         minl= g%minmaxlonlat(1,1)
         maxl= g%minmaxlonlat(1,1)
         DO i=1,SIZE(test)
            IF (test(i) > 180._dp) THEN
               test(i) = test(i) -360._dp
            END IF
            minl = MIN(test(i),minl)
            maxl = MAX(test(i),maxl)
         END DO
         diff2 = maxl-minl
         ! There must be a gap between max and 180._dp
         IF (180._dp - maxl > 15._dp) THEN
            ! FLIP RANGE
            ! FLIP LONGITUDE TO [-180:180]
            DO ix = 1, SIZE( g%ulonm%dat%vd)
               IF (g%ulonm%dat%vd(ix) > 180._dp) THEN
                  g%ulonm%dat%vd(ix) = g%ulonm%dat%vd(ix) - 360._dp
               END IF
            END DO
            DO ix = 1, SIZE( g%uloni%dat%vd)
               IF (g%uloni%dat%vd(ix) > 180._dp) THEN
                  g%uloni%dat%vd(ix) = g%uloni%dat%vd(ix) - 360._dp
               END IF
            END DO
            g%minmaxlonlat(1,1) = minl
            g%minmaxlonlat(1,2) = maxl
            write (0,*) 'FLIPPING1 ',  g%minmaxlonlat(1,:)
         ELSE
            g%minmaxlonlat(1,:) = RGEMPTY
            write (0,*) 'NOT FLIPPING1', minl,maxl,  g%minmaxlonlat(1,:)
         END IF
         DEALLOCATE(test)

      ELSE
         ! LARGE (> 180 deg) LONGITUDINAL GRID ON PE CROSSING [360:0] border
         ! i.e. gap smaller than 180 deg and thus no flipping possible
         ! => use full input grid
         ! (!! THIS IS NOT NICE BECAUSE ITS INCREASING THE REQUIRED
         !     COMPUTING TIME !)
         !
!!$       IF  (g%minmaxlonlat(1,2)-g%minmaxlonlat(1,1) > 180._dp) THEN
!!$          g%minmaxlonlat(1,:) = RGEMPTY
!!$          write (0,*) 'GRID_BI LON corr: ', g%minmaxlonlat(1,:)
!!$       END IF
         !      write (0,*) 'NOT FLIPPING ',  g%minmaxlonlat(1,:)
      END IF

      NULLIFY(gdat)
      NULLIFY(ldat)

      DEALLOCATE(dum2dloc)
      NULLIFY(dum2dloc)
      DEALLOCATE(dummy2d)
      NULLIFY(dummy2d)
      CALL INIT_NCVAR(loclon)

      ! -------------------------------------
      ! Unstructured LATITUDE corners (INTERFACES) ...
      g%ulati%name  = 'unstruc_corners_lat'
      g%ulati%id    = NULL_VARID
      g%ulati%xtype = NF90_DOUBLE
      ! ... dimensions
      g%ulati%ndims = 3
      ALLOCATE(g%ulati%dim(g%ulati%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%ulati%dim(1)%name = 'unstruc_lon'
      g%ulati%dim(1)%id    = NULL_DIMID
      g%ulati%dim(1)%len   = g%ulatm%dim(1)%len
      g%ulati%dim(1)%fuid  = .false.
      g%ulati%dim(1)%varid = NULL_VARID
      g%ulati%dim(2)%name = 'unstruc_lat'
      g%ulati%dim(2)%id    = NULL_DIMID
      g%ulati%dim(2)%len   = g%ulatm%dim(2)%len
      g%ulati%dim(2)%fuid  = .false.
      g%ulati%dim(2)%varid = NULL_VARID
      g%ulati%dim(3)%name  = 'num_corners'
      g%ulati%dim(3)%id    = NULL_DIMID
      g%ulati%dim(3)%len   = g%corners ! or loclat%dim(1)%len
      g%ulati%dim(3)%fuid  = .false.
      g%ulati%dim(3)%varid = NULL_VARID

      ! ... ATTRIBUTES
      g%ulati%natts = 1
      ALLOCATE(g%ulati%att(g%ulati%natts))
      DO ix = 1, loclat%natts
         IF (TRIM(loclat%att(ix)%name) == 'units') THEN
            CALL COPY_NCATT(g%ulati%att(1),loclat%att(ix))
            EXIT
         END IF
      END DO

      ! ... data
      CALL INIT_NARRAY(g%ulati%dat, g%ulati%ndims &
           , (/g%ulati%dim(1)%len,g%ulati%dim(2)%len,g%ulati%dim(3)%len/) &
           ,VTYPE_DOUBLE)

      CALL INIT_NCVAR(loclat)
      CALL IMPORT_NCVAR(loclat,varname='grid_corner_lat', file=ncdata)
      CALL DOUBLE_NARRAY(loclat%dat)
      ! convert from global 1D to  global (corners,ncol)
      CALL DEALINE_ARRAY(status,loclat%dim(1)%len, loclat%dim(2)%len &
           , loclat%dat%vd, dummy2d)

      ! convert from  global (corners,ncol) to local (nproma, ngpblks, corners)
      ALLOCATE(dum2dloc(nproma,ngpblks,5,1))
      DO ix = 1,SIZE(dummy2d,1)
         gdat => dummy2d(ix,:)
         ldat => dum2dloc(:,:,ix,1)
         CALL SCATTER_GP(gdat,ldat,dcg)
      END DO

      ! FLAG empty points
      DO jx = 1, ngpblks
         IF (npromz(jx) /= nproma) THEN
            DO ix = npromz(jx)+1, nproma
               dum2dloc(ix,jx,:,:) = dum2dloc(npromz(jx),jx,:,:)
            END DO
         END IF
      END DO
      g%minmaxlonlat(2,1) = MINVAL(dum2dloc)
      g%minmaxlonlat(2,2) = MAXVAL(dum2dloc)

      IF (MAXVAL(ABS(g%minmaxlonlat(2,:))) > 89._dp) THEN
         g%minmaxlonlat(:,:) = RGEMPTY
      END IF
      ! ... assign 1d GRID Data
      dim1=g%ulati%dim(1)%len
      gsize=g%ulati%dim(1)%len *  g%ulati%dim(2)%len
      DO n = 1, g%corners
         DO ix = 1, g%ulati%dim(1)%len
            DO jx = 1, g%ulati%dim(2)%len
               g%ulati%dat%vd((jx-1)*dim1 + ix + (n-1)*gsize) = &
                    dum2dloc(ix,jx,n,1)
            END DO
         END DO
      END DO

      NULLIFY(gdat)
      NULLIFY(ldat)

      DEALLOCATE(dum2dloc)
      NULLIFY(dum2dloc)
      DEALLOCATE(dummy2d)
      NULLIFY(dummy2d)
      CALL INIT_NCVAR(loclat)
      CALL INIT_NCVAR(loclon)

      ! -------------------------------

      g%lonc = .FALSE. ! Axis is not circular
#endif

      ! HYBRID-A-COEFFICIENTS (INTERFACES) ...
      g%hyai%name  = 'hyai'
      g%hyai%id    = NULL_VARID
      g%hyai%xtype = NF90_FLOAT
      ! ... dimensions
      g%hyai%ndims = 1
      ALLOCATE(g%hyai%dim(g%hyai%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
      g%hyai%dim(1)%name  = 'ilev'
      g%hyai%dim(1)%id    = NULL_DIMID
      g%hyai%dim(1)%len   = nlev+1
      g%hyai%dim(1)%fuid  = .false.
      g%hyai%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%hyai%dat, g%hyai%ndims, (/g%hyai%dim(1)%len/) &
           ,VTYPE_REAL)
      g%hyai%dat%vr(:) = REAL(vct(1:nvclev),sp)
      ! ... attributes
      CALL ADD_NCATT(g%hyai, 'long_name'                              &
           ,vs='hybrid-A-coefficients at layer interfaces')
      CALL ADD_NCATT(g%hyai, 'units', vs='1')
      g%ranges(3,1) = 0.0_dp
      g%ranges(3,2) = 0.0_dp

      ! HYBRID-B-COEFFICIENTS (INTERFACES) ...
      g%hybi%name  = 'hybi'
      g%hybi%id    = NULL_VARID
      g%hybi%xtype = NF90_FLOAT
      ! ... dimensions
      g%hybi%ndims = 1
      ALLOCATE(g%hybi%dim(g%hybi%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,4)
      g%hybi%dim(1) = g%hyai%dim(1)
      ! ... data
      CALL INIT_NARRAY(g%hybi%dat, g%hybi%ndims, (/g%hybi%dim(1)%len/) &
           ,VTYPE_REAL)
      g%hybi%dat%vr(:) = REAL(vct(nvclev+1:2*nvclev),sp)
      ! ... attributes
      CALL ADD_NCATT(g%hybi, 'long_name'                              &
           ,vs='hybrid-B-coefficients at layer interfaces')
      CALL ADD_NCATT(g%hybi, 'units', vs='1')
      g%ranges(4,1) = 0.0_dp
      g%ranges(4,2) = 1.0_dp

      ! SURFACE PRESSURE
      g%ps%name  = 'ps'
      g%ps%id    = NULL_VARID
      g%ps%xtype = NF90_FLOAT
      ! ... dimensions
      g%ps%ndims = 2
      ALLOCATE(g%ps%dim(g%ps%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)
#ifndef HOMMESE
      g%ps%dim(1) = g%lonm%dim(1)
      g%ps%dim(2) = g%latm%dim(1)
#else
      g%ps%dim(1) = g%ulonm%dim(1)
      g%ps%dim(2) = g%ulonm%dim(2)
#endif
      ! ... data
      CALL INIT_NARRAY(g%ps%dat, g%ps%ndims                  &
           ,(/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
           ,VTYPE_REAL)
      IF (ASSOCIATED(aps)) THEN
         DO i=1, g%ps%dim(1)%len
            DO j=1, g%ps%dim(2)%len
               g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/)  &
                    ,(/i,j/)))                           &
                    !qqq
                    ! = aps(i,j)
                    = 101325.0
            END DO
         END DO
      ELSE
         DO i=1, g%ps%dim(1)%len
            DO j=1, g%ps%dim(2)%len
               g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/)  &
                    ,(/i,j/)))                           &
                    = 101325.0
            END DO
         END DO
      END IF
      ! ... attributes
      CALL ADD_NCATT(g%ps, 'long_name', vs='surface pressure')
      CALL ADD_NCATT(g%ps, 'units', vs='Pa')

      ! REFERENCE PRESSURE ...
      g%p0%name  = 'p0'
      g%p0%id    = NULL_VARID
      g%p0%xtype = NF90_FLOAT
      ! ... dimensions
      g%p0%ndims = 0
      ! ... data
      CALL INIT_NARRAY(g%p0%dat, 1, (/ 1 /), VTYPE_REAL)
      g%p0%dat%vr(1) = 1.0
      ! ... attributes
      CALL ADD_NCATT(g%p0, 'long_name', vs='reference pressure')
      CALL ADD_NCATT(g%p0, 'units', vs='Pa')

      ! TIME (MID) ...
      g%timem%name  = 'time'
      g%timem%id    = NULL_VARID
      g%timem%xtype = NF90_FLOAT
      ! ... dimensions
      g%timem%ndims = 1
      ALLOCATE(g%timem%dim(g%timem%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
      g%timem%dim(1)%name  = 'time'
      g%timem%dim(1)%id    = NULL_DIMID
      g%timem%dim(1)%len   = 1
      g%timem%dim(1)%fuid  = .true.
      g%timem%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%timem%dat, g%timem%ndims, (/g%timem%dim(1)%len/) &
           ,VTYPE_DOUBLE)
      ! TIME: SECONDS SINCE MODEL START
      CALL TIME_SPAN_D(dts   &
           , YEAR_START, MONTH_START, DAY_START &
           , HOUR_START, MINUTE_START, SECOND_START  &
           , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
      g%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds
      ! ... attributes
      CALL ADD_NCATT(g%timem, 'long_name'                              &
           ,vs='time in seconds since model start')
      WRITE(tunit, &
           '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
           YEAR_START, MONTH_START, DAY_START &
           , HOUR_START, MINUTE_START, SECOND_START
      CALL ADD_NCATT(g%timem, 'units', vs=TRIM(tunit))

      ! CALCULATE INTs from MIDs
      CALL COMPLETE_GEOHYBGRID(g)

      ! INTERFACE OK
      ok = .true.

      !    CALL PRINT_GEOHYBGRID(g, 'GRID_BI',  60)

    END SUBROUTINE DEFINE_BASEMODEL_GEOHYBGRID
    ! -------------------------------------------------------------------------
#endif
!CESM1

!ICON
#ifdef ICON
    ! -------------------------------------------------------------------------
    SUBROUTINE DEFINE_BASEMODEL_GEOHYBGRID(g, ok)

      USE messy_main_bmluse_bi,     ONLY: p_patch, nproma
      USE messy_main_timer,         ONLY: YEAR_START, MONTH_START, DAY_START &
                                        , HOUR_START, MINUTE_START &
                                        , SECOND_START             &
                                        , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
      USE messy_main_timer,         ONLY: TIME_SPAN_D
      USE messy_main_tools,         ONLY: int2str
      USE messy_main_grid_trafo,    ONLY: RGEMPTY, RGMAX, COMPLETE_GEOHYBGRID
      USE messy_main_grid_tools,    ONLY: RGTOOL_CONVERT_DAT2PREDEFVAR
      USE messy_main_grid_netcdf,   ONLY: ERRMSG, NULL_DIMID, NULL_VARID &
                                        , VTYPE_DOUBLE, VTYPE_INT        &
                                        , nf90_float, NF90_DOUBLE        &
                                        , INIT_NARRAY, ADD_NCATT         &
                                        , INIT_NCVAR
      USE messy_main_grid,          ONLY: t_geohybgrid, INIT_GEOHYBGRID
      USE messy_main_channel_mem,   ONLY: dom_curid
      USE messy_main_constants_mem, ONLY: RTD_icon
      USE mo_impl_constants,        ONLY: min_rlcell
      USE mo_loopindices,           ONLY: get_indices_c

      IMPLICIT NONE

      ! I/O
      TYPE (t_geohybgrid), INTENT(OUT) :: g
      LOGICAL            , INTENT(OUT) :: ok

      ! LOCAL
      CHARACTER(LEN=2)    :: spid                ! string dom_id
      INTEGER             :: vertind1, vertind2
      INTEGER             :: n, ix,jx, dim1, gsize
      REAL(DP)            :: dts
      INTEGER             :: status
      CHARACTER(LEN=100)  :: tunit
      INTEGER             :: icp, is, ie, jc, jb
      INTEGER             :: i_startblk, i_endblk
      REAL(DP)            :: center_lon, right_bound

      icp = dom_curid

      ! POINTER FOR SCATTER ROUTINES
      ! INIT
      CALL INIT_GEOHYBGRID(g)

      CALL int2str(spid,icp,'0')
      g%name    = 'ICONp'//spid
      g%file    = ''                              ! Filename
      g%t       = 0                               ! time step
      g%corners = 3                               ! number of corners

      ! unstructured grid
      g%lonc = .FALSE. ! Axis is not circular

      ! DEFINE MASK of defined points
      g%imask%name = 'mask'
      g%imask%id    = NULL_VARID
      g%imask%xtype = NF90_DOUBLE
      ! ... dimensions
      g%imask%ndims = 2
      ALLOCATE(g%imask%dim(g%imask%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%imask%dim(1)%name = 'unstruc_lon'
      g%imask%dim(1)%id    = NULL_DIMID
      g%imask%dim(1)%len   = nproma
      g%imask%dim(1)%fuid  = .false.
      g%imask%dim(1)%varid = NULL_VARID
      g%imask%dim(2)%name  = 'unstruc_lat'
      g%imask%dim(2)%id    = NULL_DIMID
      g%imask%dim(2)%len   = p_patch(icp)%nblks_c
      g%imask%dim(2)%fuid  = .false.
      g%imask%dim(2)%varid = NULL_VARID

      ! .. ATTRIBUTES
      CALL ADD_NCATT(g%imask, 'long_name', vs='mask of defined points')

      ! ... data
      CALL INIT_NARRAY(g%imask%dat, g%imask%ndims &
           , (/g%imask%dim(1)%len,g%imask%dim(2)%len/) &
           ,VTYPE_INT)

      g%imask%dat%vi = 1

      i_startblk = p_patch(icp)%cells%start_block(1)
      i_endblk   = p_patch(icp)%cells%end_block(min_rlcell)

      DO jc = 1, nproma
         DO jb = 1, p_patch(icp)%nblks_c
            IF (jb < i_startblk .OR. jb > i_endblk) THEN
               g%imask%dat%vi((jb-1)* nproma + jc) = 0
            ELSE
               CALL get_indices_c(p_patch(icp), jb, i_startblk,i_endblk &
                    , is, ie &
                    , 1, min_rlcell)
               IF (jc < is .OR. jc > ie) THEN
                  g%imask%dat%vi((jb-1)*nproma  + jc) = 0
               END IF
            END IF
         END DO
      END DO

      ! Unstructured LONGITUDE Centers(MID) ...
      CALL INIT_NCVAR(g%ulonm)
      CALL INIT_NCVAR(g%ulatm)

      g%ulonm%name  = 'unstruc_lon'
      g%ulonm%id    = NULL_VARID
      g%ulonm%xtype = NF90_DOUBLE
      ! ... dimensions
      g%ulonm%ndims = 2

      ALLOCATE(g%ulonm%dim(g%ulonm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%ulonm%dim(1)%name = 'unstruc_lon'
      g%ulonm%dim(1)%id    = NULL_DIMID
      g%ulonm%dim(1)%len   = nproma
      g%ulonm%dim(1)%fuid  = .false.
      g%ulonm%dim(1)%varid = NULL_VARID


      g%ulonm%dim(2)%name  = 'unstruc_lat'
      g%ulonm%dim(2)%id    = NULL_DIMID
      g%ulonm%dim(2)%len   = p_patch(icp)%nblks_c
      g%ulonm%dim(2)%fuid  = .false.
      g%ulonm%dim(2)%varid = NULL_VARID
      g%ranges(1,1) = RGEMPTY
      g%ranges(1,2) = RGEMPTY

      ! .. ATTRIBUTES
      CALL ADD_NCATT(g%ulonm, 'long_name', vs='longitude')
      CALL ADD_NCATT(g%ulonm, 'units', vs='degrees_east')

      CALL INIT_NARRAY(g%ulonm%dat, g%ulonm%ndims &
           , (/g%ulonm%dim(1)%len,g%ulonm%dim(2)%len/) &
           , VTYPE_DOUBLE)

      ! ... assign 1d GRID Data
      dim1=g%ulonm%dim(1)%len
      DO ix = 1, g%ulonm%dim(1)%len
         DO jx = 1, g%ulonm%dim(2)%len
            IF (g%imask%dat%vi((jx-1)*dim1 + ix) > 0) THEN
               g%ulonm%dat%vd((jx-1)*dim1 + ix) = &
                    p_patch(icp)%cells%center(ix,jx)%lon * RTD_icon
            ELSE
               g%ulonm%dat%vd((jx-1)*dim1 + ix) = -999._dp
            END IF
         END DO
      END DO

      ! -------------------------------------
      ! Unstructured LATITUDE Centers(MID) ...
      g%ulatm%name  = 'unstruc_lat'
      g%ulatm%id    = NULL_VARID
      g%ulatm%xtype = NF90_DOUBLE
      ! ... dimensions
      g%ulatm%ndims = 2
      ALLOCATE(g%ulatm%dim(g%ulatm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%ulatm%dim(1)%name = 'unstruc_lon'
      g%ulatm%dim(1)%id    = NULL_DIMID
      g%ulatm%dim(1)%len   = nproma
      g%ulatm%dim(1)%fuid  = .false.
      g%ulatm%dim(1)%varid = NULL_VARID

      g%ulatm%dim(2)%name  = 'unstruc_lat'
      g%ulatm%dim(2)%id    = NULL_DIMID
      g%ulatm%dim(2)%len   = p_patch(dom_curid)%nblks_c
      g%ulatm%dim(2)%fuid  = .false.
      g%ulatm%dim(2)%varid = NULL_VARID

      g%ranges(2,1) = RGEMPTY !-90.0_dp
      g%ranges(2,2) = RGEMPTY ! 90.0_dp
      ! .. ATTRIBUTES
      CALL ADD_NCATT(g%ulatm, 'long_name', vs='latitude')
      CALL ADD_NCATT(g%ulatm, 'units', vs='degrees_north')

      ! ... data
      CALL INIT_NARRAY(g%ulatm%dat, g%ulatm%ndims &
           , (/g%ulatm%dim(1)%len,g%ulatm%dim(2)%len/) &
           ,VTYPE_DOUBLE)

      ! ... assign 1d GRID Dat
      dim1=g%ulatm%dim(1)%len
      DO ix = 1, g%ulatm%dim(1)%len
         DO jx = 1, g%ulatm%dim(2)%len
            IF (g%imask%dat%vi((jx-1)*dim1 + ix) > 0) THEN
               g%ulatm%dat%vd((jx-1)*dim1 + ix) = &
                    p_patch(icp)%cells%center(ix,jx)%lat * RTD_icon
            ELSE
               g%ulatm%dat%vd((jx-1)*dim1 + ix) = -999.0_dp
            END IF
         END DO
      END DO
      ! -------------------------------------

      ! -------------------------------------
      CALL INIT_NCVAR(g%uloni)
      CALL INIT_NCVAR(g%ulati)
      ! -------------------------------------

      ! -------------------------------------
      ! Unstructured LONGITUDE corners (INTERFACES) ...
      g%uloni%name  = 'unstruc_corners_lon'
      g%uloni%id    = NULL_VARID
      g%uloni%xtype = NF90_DOUBLE
      ! ... dimensions
      g%uloni%ndims = 3
      ALLOCATE(g%uloni%dim(g%uloni%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%uloni%dim(1)%name = 'unstruc_lon'
      g%uloni%dim(1)%id    = NULL_DIMID
      g%uloni%dim(1)%len   = g%ulonm%dim(1)%len
      g%uloni%dim(1)%fuid  = .false.
      g%uloni%dim(1)%varid = NULL_VARID
      g%uloni%dim(2)%name = 'unstruc_lon'
      g%uloni%dim(2)%id    = NULL_DIMID
      g%uloni%dim(2)%len   = g%ulonm%dim(2)%len
      g%uloni%dim(2)%fuid  = .false.
      g%uloni%dim(2)%varid = NULL_VARID

      g%uloni%dim(3)%name  = 'num_corners'
      g%uloni%dim(3)%id    = NULL_DIMID
      g%uloni%dim(3)%len   = g%corners ! or loclon%dim(1)%len
      g%uloni%dim(3)%fuid  = .false.
      g%uloni%dim(3)%varid = NULL_VARID

      ! COPY ATTRIBUTES
      CALL ADD_NCATT(g%uloni, 'long_name', vs='corner longitude')
      CALL ADD_NCATT(g%uloni, 'units',     vs='degrees_east')

      ! ... data
      CALL INIT_NARRAY(g%uloni%dat, g%uloni%ndims &
           , (/g%uloni%dim(1)%len,g%uloni%dim(2)%len,g%uloni%dim(3)%len/) &
           ,VTYPE_DOUBLE)

      ! ... assign 1d GRID Data
      dim1=g%uloni%dim(1)%len
      gsize=g%uloni%dim(1)%len *  g%uloni%dim(2)%len
      DO n = 1, g%corners
         DO ix = 1, g%uloni%dim(1)%len
            DO jx = 1, g%uloni%dim(2)%len
               IF (jx < i_startblk .OR. jx > i_endblk) THEN
                  g%uloni%dat%vd((jx-1)*dim1 + ix + (n-1)*gsize) = -999.0_dp
               ELSE
                  CALL get_indices_c(p_patch(icp),jx, i_startblk &
                       , i_endblk, is, ie &
                       , 1, min_rlcell)
                  IF (ix < is .OR. ix > ie) THEN
                     g%uloni%dat%vd((jx-1)*dim1 + ix + (n-1)*gsize) = -999.0_dp
                  ELSE
                     vertind1 = p_patch(icp)%cells%vertex_idx(ix,jx,n)
                     vertind2 = p_patch(icp)%cells%vertex_blk(ix,jx,n)

                     g%uloni%dat%vd((jx-1)*dim1 + ix + (n-1)*gsize) = &
                          p_patch(icp)%verts%vertex(vertind1,vertind2)%lon &
                          * RTD_icon
                  END IF
               END IF
            END DO
         END DO
      END DO

      g%minmaxlonlat(1,1) = MINVAL(g%uloni%dat%vd)
      g%minmaxlonlat(1,2) = MAXVAL(g%uloni%dat%vd)

      ! Test if the domain spans over the grid's boundary and set the
      ! correct MIN/MAX(lon) values.
      ! Crucial test, if domain spans over half the globe. This is unliekely,
      ! so wie may deal with a modulo domain, i.e. wrapping around
      ! at 0° or 180°
      IF (g%minmaxlonlat(1,2) - g%minmaxlonlat(1,1) > 180._dp) THEN
         ! initialise center and right boundary of the grid (0° - 360°)
         center_lon  = 180._dp
         right_bound = 360._dp
         ! test for domain on (-180° - 180°)
         IF (ANY(g%uloni%dat%vd < 0._dp)) THEN
            center_lon = 0.
            right_bound = 180._dp
         END IF
         ! wrap the lower half of the grid to higher end and find new min/max
         WHERE (g%uloni%dat%vd < center_lon) &
              g%uloni%dat%vd = g%uloni%dat%vd + 360._dp
         g%minmaxlonlat(1,1) = MINVAL(g%uloni%dat%vd)
         g%minmaxlonlat(1,2) = MAXVAL(g%uloni%dat%vd)
         ! shift back and move found maximum in original interval
         WHERE (g%uloni%dat%vd > right_bound) &
              g%uloni%dat%vd = g%uloni%dat%vd - 360._dp
         g%minmaxlonlat(1,2) = g%minmaxlonlat(1,2) - 360._dp
         ! test if the domain was really wrapped around,
         ! or if it is in original grid interval
         IF (ANY((g%uloni%dat%vd > g%minmaxlonlat(1,2))              &
              &    .AND. (g%uloni%dat%vd < g%minmaxlonlat(1,1)))) THEN
            g%minmaxlonlat(1,1) = MINVAL(g%uloni%dat%vd)
            g%minmaxlonlat(1,2) = MAXVAL(g%uloni%dat%vd)
         END IF
         ! at this point, MIN(lon) < MAX(lon) for "normal" domain,
         ! MIN(lon) > MAX(lon) for domain spanning over the grid's boundary in x
      END IF

      ! -------------------------------------
      ! Unstructured LATITUDE corners (INTERFACES) ...
      g%ulati%name  = 'unstruc_corners_lat'
      g%ulati%id    = NULL_VARID
      g%ulati%xtype = NF90_DOUBLE
      ! ... dimensions
      g%ulati%ndims = 3
      ALLOCATE(g%ulati%dim(g%ulati%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%ulati%dim(1)%name = 'unstruc_lon'
      g%ulati%dim(1)%id    = NULL_DIMID
      g%ulati%dim(1)%len   = g%ulatm%dim(1)%len
      g%ulati%dim(1)%fuid  = .false.
      g%ulati%dim(1)%varid = NULL_VARID
      g%ulati%dim(2)%name = 'unstruc_lat'
      g%ulati%dim(2)%id    = NULL_DIMID
      g%ulati%dim(2)%len   = g%ulatm%dim(2)%len
      g%ulati%dim(2)%fuid  = .false.
      g%ulati%dim(2)%varid = NULL_VARID
      g%ulati%dim(3)%name  = 'num_corners'
      g%ulati%dim(3)%id    = NULL_DIMID
      g%ulati%dim(3)%len   = g%corners ! or loclat%dim(1)%len
      g%ulati%dim(3)%fuid  = .false.
      g%ulati%dim(3)%varid = NULL_VARID

      ! ... ATTRIBUTES
      CALL ADD_NCATT(g%ulati, 'long_name', vs='corner latitude')
      CALL ADD_NCATT(g%ulati, 'units',     vs='degrees_north')

      ! ... data
      CALL INIT_NARRAY(g%ulati%dat, g%ulati%ndims &
           , (/g%ulati%dim(1)%len,g%ulati%dim(2)%len,g%ulati%dim(3)%len/) &
           ,VTYPE_DOUBLE)

      ! ... assign 1d GRID Data
      dim1=g%ulati%dim(1)%len
      gsize=g%ulati%dim(1)%len *  g%ulati%dim(2)%len
      DO n = 1, g%corners
         DO ix = 1, g%ulati%dim(1)%len
            DO jx = 1, g%ulati%dim(2)%len
               IF (jx < i_startblk .OR. jx > i_endblk) THEN
                  g%ulati%dat%vd((jx-1)*dim1 + ix + (n-1)*gsize) = -999.0_dp
               ELSE
                  CALL get_indices_c(p_patch(icp),jx, i_startblk &
                       , i_endblk, is, ie &
                       , 1, min_rlcell)
                  IF (ix < is .OR. ix > ie) THEN
                     g%ulati%dat%vd((jx-1)*dim1 + ix + (n-1)*gsize) = -999.0_dp
                  ELSE
                     vertind1 = p_patch(icp)%cells%vertex_idx(ix,jx,n)
                     vertind2 = p_patch(icp)%cells%vertex_blk(ix,jx,n)
                     g%ulati%dat%vd((jx-1)*dim1 + ix + (n-1)*gsize) = &
                          p_patch(icp)%verts%vertex(vertind1,vertind2)%lat &
                          * RTD_icon
                  END IF
               END IF
            END DO
         END DO
      END DO

      g%minmaxlonlat(2,1) = MINVAL(g%ulati%dat%vd &
           , mask=g%ulati%dat%vd > -999.0_dp)
      g%minmaxlonlat(2,2) = MAXVAL(g%ulati%dat%vd)

      IF (g%minmaxlonlat(1,2) == RGEMPTY) then
         g%minmaxlonlat(2,:) = RGEMPTY
      END IF

      ! -------------------------------

      g%lonc = .FALSE. ! Axis is not circular

      ! 3D pressure field at vertical interface levels
      g%pressi%name  = 'pressi'
      g%pressi%id    = NULL_VARID
      g%pressi%xtype = NF90_DOUBLE
      ! ... dimensions
      g%pressi%ndims = 3
      ALLOCATE(g%pressi%dim(g%pressi%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
      g%pressi%dim(1) = g%ulonm%dim(1)
      g%pressi%dim(2) = g%ulonm%dim(2)

      g%pressi%dim(3)%name  = 'ilev'
      g%pressi%dim(3)%id    = NULL_DIMID
      g%pressi%dim(3)%len   = p_patch(icp)%nlev+1
      g%pressi%dim(3)%fuid  = .false.
      g%pressi%dim(3)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%pressi%dat, g%pressi%ndims &
           , (/g%pressi%dim(1)%len,g%pressi%dim(2)%len,g%pressi%dim(3)%len/) &
           , VTYPE_DOUBLE)
      ! ... attributes
      CALL ADD_NCATT(g%pressi, 'long_name'                              &
           ,vs='3d pressure at layer interfaces')
      CALL ADD_NCATT(g%pressi, 'units', vs='Pa')
      g%ranges(3,1) = 0.0_dp
      g%ranges(3,2) = RGMAX

      ! 3D pressure field at vertical mid levels
      g%pressm%name  = 'pressm'
      g%pressm%id    = NULL_VARID
      g%pressm%xtype = NF90_DOUBLE
      ! ... dimensions
      g%pressm%ndims = 3
      ALLOCATE(g%pressm%dim(g%pressm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
      g%pressm%dim(1) = g%ulonm%dim(1)
      g%pressm%dim(2) = g%ulonm%dim(2)

      g%pressm%dim(3)%name  = 'ilev'
      g%pressm%dim(3)%id    = NULL_DIMID
      g%pressm%dim(3)%len   = p_patch(icp)%nlev
      g%pressm%dim(3)%fuid  = .false.
      g%pressm%dim(3)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%pressm%dat, g%pressm%ndims &
           , (/g%pressm%dim(1)%len,g%pressm%dim(2)%len,g%pressm%dim(3)%len/) &
           , VTYPE_DOUBLE)
      ! ... attributes
      CALL ADD_NCATT(g%pressm, 'long_name'                              &
           ,vs='3d pressure at layer interfaces')
      CALL ADD_NCATT(g%pressm, 'units', vs='Pa')

      ! TIME (MID) ...
      g%timem%name  = 'time'
      g%timem%id    = NULL_VARID
      g%timem%xtype = NF90_FLOAT
      ! ... dimensions
      g%timem%ndims = 1
      ALLOCATE(g%timem%dim(g%timem%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
      g%timem%dim(1)%name  = 'time'
      g%timem%dim(1)%id    = NULL_DIMID
      g%timem%dim(1)%len   = 1
      g%timem%dim(1)%fuid  = .true.
      g%timem%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%timem%dat, g%timem%ndims, (/g%timem%dim(1)%len/) &
           ,VTYPE_DOUBLE)
      ! TIME: SECONDS SINCE MODEL START
      CALL TIME_SPAN_D(dts   &
           , YEAR_START, MONTH_START, DAY_START &
           , HOUR_START, MINUTE_START, SECOND_START  &
           , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
      g%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds
      ! ... attributes
      CALL ADD_NCATT(g%timem, 'long_name'                              &
           ,vs='time in seconds since model start')
      WRITE(tunit, &
           '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
           YEAR_START, MONTH_START, DAY_START &
           , HOUR_START, MINUTE_START, SECOND_START
      CALL ADD_NCATT(g%timem, 'units', vs=TRIM(tunit))

      ! CALCULATE INTs from MIDs
      CALL COMPLETE_GEOHYBGRID(g)

      ! INTERFACE OK
      ok = .true.

    END SUBROUTINE DEFINE_BASEMODEL_GEOHYBGRID
    ! -------------------------------------------------------------------------
#endif
!ICON

#if defined(BLANK)
    ! -------------------------------------------------------------------------
    SUBROUTINE DEFINE_BASEMODEL_GEOHYBGRID(g, ok)

      USE messy_main_grid_def_mem_bi, ONLY: nlon, nlev, nlat
      USE messy_main_grid_def_bi,     ONLY: longitude, latitude &
#if defined(MBM_QBO)
                                          , hyai=> ceta
#elif defined(MBM_RAD)
                                          , hyai, hybi, hyam &
                                          , nvclev, vct, lv_echam
#else
                                          , hyai
#endif
      USE messy_main_timer,        ONLY: YEAR_START, MONTH_START, DAY_START &
                                       , HOUR_START, MINUTE_START           &
                                       , SECOND_START                       &
                                       , YEAR, MONTH, DAY                   &
                                       , HOUR, MINUTE, SECOND
      USE messy_main_constants_mem, ONLY: sp
      USE messy_main_timer,         ONLY: TIME_SPAN_D
      USE messy_main_grid_trafo,    ONLY: RGEMPTY, COMPLETE_GEOHYBGRID
      USE messy_main_grid_netcdf,   ONLY: ERRMSG                              &
                                        , NULL_DIMID, NULL_VARID              &
                                        , VTYPE_REAL, VTYPE_DOUBLE, VTYPE_INT &
                                        , NF90_FLOAT, NF90_DOUBLE             &
                                        , POSITION, INIT_NARRAY, ADD_NCATT
      USE messy_main_grid,          ONLY: t_geohybgrid, INIT_GEOHYBGRID

      IMPLICIT NONE

      ! I/O
      TYPE(t_geohybgrid), INTENT(OUT) :: g
      LOGICAL           , INTENT(OUT) :: ok

      ! LOCAL
      REAL(DP)                  :: dts
      INTEGER                   :: i,j
      INTEGER                   :: status
      CHARACTER(LEN=100)        :: tunit
      REAL(dp)                  :: dx

      ! INIT
      CALL INIT_GEOHYBGRID(g)

      g%name    = 'BLANK'
      g%file    = ''                           ! Filename
      g%t       = 0                            ! time step
      g%corners = 4                            ! number of corners

#if defined(MBM_RAD)
      g%lonc = .TRUE.  ! Axis is circular
      g%name = 'MBM_RAD'
#else
      g%lonc = .FALSE. ! Axis is not circular
#endif
      ! Curvilinear LONGITUDE (MID) ...
      g%lonm%name  = 'lon'
      g%lonm%id    = NULL_VARID
      g%lonm%xtype = NF90_DOUBLE
      ! ... dimensions
      g%lonm%ndims = 1
      ALLOCATE(g%lonm%dim(g%lonm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%lonm%dim(1)%name  = 'lon'
      g%lonm%dim(1)%id    = NULL_DIMID
      g%lonm%dim(1)%len   = nlon
      g%lonm%dim(1)%fuid  = .false.
      g%lonm%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%lonm%dat, g%lonm%ndims &
           , (/g%lonm%dim(1)%len/),VTYPE_DOUBLE)

      DO i=1, nlon
         g%lonm%dat%vd(i) = longitude(i)
      END DO

      ! LATITUDE (MID) ...
      g%latm%name  = 'lat'
      g%latm%id    = NULL_VARID
      g%latm%xtype = NF90_DOUBLE
      ! ... dimensions
      g%latm%ndims = 1
      ALLOCATE(g%latm%dim(g%latm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
      g%latm%dim(1)%name  = 'lat'
      g%latm%dim(1)%id    = NULL_DIMID
      g%latm%dim(1)%len   = nlat
      g%latm%dim(1)%fuid  = .false.
      g%latm%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%latm%dat, g%latm%ndims &
           , (/g%latm%dim(1)%len/),VTYPE_DOUBLE)

      DO i=1, nlat
         g%latm%dat%vd(i) = latitude(i)
      END DO
      ! ... attributes
#if defined(MBM_RAD)
      CALL ADD_NCATT(g%lonm, 'long_name', vs='longitude')
      CALL ADD_NCATT(g%lonm, 'units', vs='degrees_east')

      CALL ADD_NCATT(g%latm, 'long_name', vs='latitude')
      CALL ADD_NCATT(g%latm, 'units', vs='degrees_north')
#else
      CALL ADD_NCATT(g%lonm, 'long_name', vs='mid longitude')
      CALL ADD_NCATT(g%lonm, 'units',     vs='degrees_east')

      CALL ADD_NCATT(g%latm, 'long_name', vs='mid latitude')
      CALL ADD_NCATT(g%latm, 'units', vs='degrees_north')
#endif
      ! ----------------------------------------------------------
      ! define interfaces:
      ! Curvilinear LONGITUDE (INTERFACES) ...
      g%loni%name  = 'lon_I'
      g%loni%id    = NULL_VARID
      g%loni%xtype = NF90_DOUBLE
      ! ... dimensions
      g%loni%ndims = 1
      ALLOCATE(g%loni%dim(g%loni%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%loni%dim(1)%name  = 'lon_I'
      g%loni%dim(1)%id    = NULL_DIMID
      g%loni%dim(1)%len   = nlon+1
      g%loni%dim(1)%fuid  = .false.
      g%loni%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%loni%dat, g%loni%ndims &
           , (/g%loni%dim(1)%len/),VTYPE_DOUBLE)

      IF (SIZE(longitude) > 1) THEN
         ! This MIN(SIZE ...) is a NAG compiler work around
         dx =  (longitude(MIN(SIZE(longitude),2)) - longitude(1))/2._dp
         DO i=1, nlon
            g%loni%dat%vd(i) = longitude(i) - dx
         END DO
         g%loni%dat%vd(nlon+1) = longitude(nlon) + dx
      ELSE
         IF (NINT(longitude(1)) == 0) THEN
#if defined(MBM_RAD)
            g%loni%dat%vd(1) = longitude(1)
            g%loni%dat%vd(2) = 360._dp
#else
            g%loni%dat%vd(1) = -180._dp
            g%loni%dat%vd(2) = 180._dp
#endif
         ELSE
            ! JUST decide for an artificial value
            g%loni%dat%vd(1) = longitude(1)-0.5_dp
            g%loni%dat%vd(2) = longitude(1)+0.5_dp
         END IF
      END IF

      ! LATITUDE (MID) ...
      g%lati%name  = 'lat_I'
      g%lati%id    = NULL_VARID
      g%lati%xtype = NF90_DOUBLE

      ! ... dimensions
      g%lati%ndims = 1
      ALLOCATE(g%lati%dim(g%lati%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
      g%lati%dim(1)%name  = 'lon_I'
      g%lati%dim(1)%id    = NULL_DIMID
      g%lati%dim(1)%len   = nlat+1
      g%lati%dim(1)%fuid  = .false.
      g%lati%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%lati%dat, g%lati%ndims &
           , (/g%lati%dim(1)%len/),VTYPE_DOUBLE)


      IF (SIZE(latitude) > 1) THEN
         ! This MIN(SIZE ...) is a NAG compiler work around
         dx =  (latitude(MIN(SIZE(latitude),2)) - latitude(1))/2.
         DO i=1, nlat
            g%lati%dat%vd(i) = latitude(i) - dx
         END DO
         g%lati%dat%vd(nlat+1) = latitude(nlat) + dx
      ELSE
         IF (NINT(latitude(1)) == 0) THEN
            g%lati%dat%vd(1) = -80._dp
            g%lati%dat%vd(2) = 80._dp
         ELSE
            ! JUST decide for an artificial value
            g%lati%dat%vd(1) = latitude(1)-0.5_dp
            g%lati%dat%vd(2) = latitude(1)+0.5_dp
         END IF
      END IF


      CALL ADD_NCATT(g%cloni, 'long_name', vs='interface longitude')
      CALL ADD_NCATT(g%cloni, 'units', vs='degrees_east')

      ! ... attributes
      CALL ADD_NCATT(g%clati, 'long_name', vs='interface latitude')
      CALL ADD_NCATT(g%clati, 'units', vs='degrees_north')

      g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
      g%ranges(1,2) = RGEMPTY
#if defined(MBM_RAD)
      g%ranges(2,1) = -90.0_dp
      g%ranges(2,2) =  90.0_dp
#else
      g%ranges(2,1) = RGEMPTY !-90.0_sp
      g%ranges(2,2) = RGEMPTY
#endif

      g%minmaxlonlat(1,1) = MINVAL(g%loni%dat%vd)
      g%minmaxlonlat(1,2) = MAXVAL(g%loni%dat%vd)

      ! DETERMIN MIN / MAX  LATITUDE
      g%minmaxlonlat(2,1) = MINVAL(g%lati%dat%vd)
      g%minmaxlonlat(2,2) = MAXVAL(g%lati%dat%vd)
      !*************************************************************************
      ! **********************************************************************
      write (*,*) 'DEFINING BASEGRID as PRESSURE GRID'

      ! HYBRID-A-COEFFICIENTS (INTERFACES) ...
      g%hyai%name  = 'hyai'
      g%hyai%id    = NULL_VARID
      g%hyai%xtype = NF90_DOUBLE
      ! ... dimensions
      g%hyai%ndims = 1
      ALLOCATE(g%hyai%dim(g%hyai%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
      g%hyai%dim(1)%name  = 'ilev'
      g%hyai%dim(1)%id    = NULL_DIMID
      g%hyai%dim(1)%len   = SIZE(hyai)
      g%hyai%dim(1)%fuid  = .false.
      g%hyai%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%hyai%dat, g%hyai%ndims, (/g%hyai%dim(1)%len/) &
           ,VTYPE_DOUBLE)
      g%hyai%dat%vd(:) = hyai(:)
      ! ... attributes
      CALL ADD_NCATT(g%hyai, 'long_name'                              &
           ,vs='hybrid-A-coefficients at layer interfaces')
      CALL ADD_NCATT(g%hyai, 'units', vs='1')
      g%ranges(3,1) = 0.0_dp
      g%ranges(3,2) = 0.0_dp
      ! ********************************************************************
      ! HYBRID-B-COEFFICIENTS (INTERFACES) ...
      g%hybi%name  = 'hybi'
      g%hybi%id    = NULL_VARID
      g%hybi%xtype = NF90_DOUBLE
      ! ... dimensions
      g%hybi%ndims = 1
      ALLOCATE(g%hybi%dim(g%hybi%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,4)
      g%hybi%dim(1) = g%hyai%dim(1)
      ! ... data
      CALL INIT_NARRAY(g%hybi%dat, g%hybi%ndims, (/g%hybi%dim(1)%len/) &
           ,VTYPE_DOUBLE)
#if defined(MBM_RAD)
      IF (lv_echam) THEN
         g%hybi%dat%vd(:) = hybi(:)
      ELSE
         g%hybi%dat%vd(:) = 0._dp
      END IF
#else
      g%hybi%dat%vd(:) = 0._dp
#endif
      ! ... attributes
      CALL ADD_NCATT(g%hybi, 'long_name'                              &
           ,vs='hybrid-B-coefficients at layer interfaces')
      CALL ADD_NCATT(g%hybi, 'units', vs='1')
      g%ranges(4,1) = 0.0_dp
      g%ranges(4,2) = 1.0_dp

#if defined(MBM_RAD)
      IF (.NOT.lv_echam) THEN
      ! HYBRID-A-COEFFICIENTS ...
      g%hyam%name  = 'hyam'
      g%hyam%id    = NULL_VARID
      g%hyam%xtype = NF90_DOUBLE
      ! ... dimensions
      g%hyam%ndims = 1
      ALLOCATE(g%hyam%dim(g%hyam%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
      g%hyam%dim(1)%name  = 'lev'
      g%hyam%dim(1)%id    = NULL_DIMID
      g%hyam%dim(1)%len   = nlev
      g%hyam%dim(1)%fuid  = .false.
      g%hyam%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%hyam%dat, g%hyam%ndims, (/g%hyam%dim(1)%len/) &
           ,VTYPE_DOUBLE)
      g%hyam%dat%vd(:) = hyam(:)
      ! ... attributes
      CALL ADD_NCATT(g%hyam, 'long_name'                              &
           ,vs='hybrid-A-coefficients at mid layer')
      CALL ADD_NCATT(g%hyam, 'units', vs='Pa')

      ! HYBRID-B COEFFICIENTS ...
      g%hybm%name  = 'hybm'
      g%hybm%id    = NULL_VARID
      g%hybm%xtype = NF90_DOUBLE
      ! ... dimensions
      g%hybm%ndims = 1
      ALLOCATE(g%hybm%dim(g%hybm%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,4)
      g%hybm%dim(1) = g%hyam%dim(1)
      ! ... data
      CALL INIT_NARRAY(g%hybm%dat, g%hybm%ndims, (/g%hybm%dim(1)%len/) &
           ,VTYPE_DOUBLE)
      g%hybm%dat%vd(:) = 0._dp

      ! ... attributes
      CALL ADD_NCATT(g%hybm, 'long_name'                              &
           ,vs='hybrid-B-coefficients at mid layer')
      CALL ADD_NCATT(g%hybm, 'units', vs='1')
      END IF
#endif

      ! ***********************************************************************
      ! SURFACE PRESSURE
      g%ps%name  = 'ps'
      g%ps%id    = NULL_VARID
      g%ps%xtype = NF90_FLOAT
      ! ... dimensions
      g%ps%ndims = 2
      ALLOCATE(g%ps%dim(g%ps%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)

      g%ps%dim(1) = g%lonm%dim(1)
      g%ps%dim(2) = g%latm%dim(1)

      ! ... data
      CALL INIT_NARRAY(g%ps%dat, g%ps%ndims                  &
           ,(/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
           ,VTYPE_REAL)
      DO i=1, g%ps%dim(1)%len
         DO j=1, g%ps%dim(2)%len
            g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/)  &
                 ,(/i,j/)))                           &
                 = 101325.0_sp
         END DO
      END DO
      ! ... attributes
      CALL ADD_NCATT(g%ps, 'long_name', vs='surface pressure')
      CALL ADD_NCATT(g%ps, 'units', vs='Pa')

      ! ***********************************************************************
      ! REFERENCE PRESSURE ...
      g%p0%name  = 'p0'
      g%p0%id    = NULL_VARID
      g%p0%xtype = NF90_FLOAT
      ! ... dimensions
      g%p0%ndims = 0
      ! ... data
      CALL INIT_NARRAY(g%p0%dat, 1, (/ 1 /), VTYPE_REAL)
      g%p0%dat%vr(1) = 1.0_sp
      ! ... attributes
      CALL ADD_NCATT(g%p0, 'long_name', vs='reference pressure')
      CALL ADD_NCATT(g%p0, 'units', vs='Pa')

      ! ***********************************************************************
      ! DEFINE MASK of defined points
      g%imask%name = 'mask'
      g%imask%id    = NULL_VARID
      g%imask%xtype = NF90_DOUBLE
      ! ... dimensions
      g%imask%ndims = 2
      ALLOCATE(g%imask%dim(g%imask%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
      g%imask%dim(1)%name = 'lon'
      g%imask%dim(1)%id    = NULL_DIMID
      g%imask%dim(1)%len   = nlon
      g%imask%dim(1)%fuid  = .false.
      g%imask%dim(1)%varid = NULL_VARID
      g%imask%dim(2)%name  = 'lat'
      g%imask%dim(2)%id    = NULL_DIMID
      g%imask%dim(2)%len   = nlat
      g%imask%dim(2)%fuid  = .false.
      g%imask%dim(2)%varid = NULL_VARID

      ! .. ATTRIBUTES
      CALL ADD_NCATT(g%imask, 'long_name', vs='mask of defined points')

      ! ... data
      CALL INIT_NARRAY(g%imask%dat, g%imask%ndims &
           , (/g%imask%dim(1)%len,g%imask%dim(2)%len/) &
           ,VTYPE_INT)

      g%imask%dat%vi(:) = 1

      ! ***********************************************************************
      ! TIME (MID) ...
      g%timem%name  = 'time'
      g%timem%id    = NULL_VARID
      g%timem%xtype = NF90_DOUBLE
      ! ... dimensions
      g%timem%ndims = 1
      ALLOCATE(g%timem%dim(g%timem%ndims), STAT=status)
      CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
      g%timem%dim(1)%name  = 'time'
      g%timem%dim(1)%id    = NULL_DIMID
      g%timem%dim(1)%len   = 1
      g%timem%dim(1)%fuid  = .true.
      g%timem%dim(1)%varid = NULL_VARID
      ! ... data
      CALL INIT_NARRAY(g%timem%dat, g%timem%ndims, (/g%timem%dim(1)%len/) &
           ,VTYPE_DOUBLE)
      ! TIME: SECONDS SINCE MODEL START
      CALL TIME_SPAN_D(dts   &
           , YEAR_START, MONTH_START, DAY_START &
           , HOUR_START, MINUTE_START, SECOND_START  &
           , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
      g%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds

      ! ... attributes
      CALL ADD_NCATT(g%timem, 'long_name'                              &
           ,vs='time in seconds since model start')
      WRITE(tunit, &
           '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
           YEAR_START, MONTH_START, DAY_START &
           , HOUR_START, MINUTE_START, SECOND_START

      CALL ADD_NCATT(g%timem, 'units', vs=TRIM(tunit))

      ! CALCULATE INTs from MIDs
      CALL COMPLETE_GEOHYBGRID(g)

      ! INTERFACE OK
      ok = .true.

    END SUBROUTINE DEFINE_BASEMODEL_GEOHYBGRID
    ! -------------------------------------------------------------------------
#endif

#if defined(MBM_CLAMS) || defined(VERTICO)
    SUBROUTINE DEFINE_BASEMODEL_GEOHYBGRID(g, ok)
      ! DUMMY ROUTINE: ATTENTION NO BASEGRID DEFINED
      !  => IMPORT NOT USABLE!

      USE messy_main_grid,          ONLY: t_geohybgrid, INIT_GEOHYBGRID
      ! I/O
      TYPE(t_geohybgrid), INTENT(OUT) :: g
      LOGICAL           , INTENT(OUT) :: ok

      CALL INIT_GEOHYBGRID(g)
      OK = .TRUE.

      write (*,*) "WARNING: no GEOHYBGRID DEFINED: IMPORT NOT USABLE!"

    END SUBROUTINE DEFINE_BASEMODEL_GEOHYBGRID
#endif

  END SUBROUTINE main_grid_init_memory
  ! ---------------------------------------------------------------------------


  ! ---------------------------------------------------------------------------
  SUBROUTINE main_grid_read_restart

    USE messy_main_grid_def_bi, ONLY: main_grid_def_read_restart

    IMPLICIT NONE

     CALL main_grid_def_read_restart

  END SUBROUTINE main_grid_read_restart
  ! ---------------------------------------------------------------------------

#ifdef COSMO
  ! ---------------------------------------------------------------------------
  SUBROUTINE main_grid_init_coupling

    USE messy_main_grid_def_bi, ONLY: main_grid_def_init_coupling

    IMPLICIT NONE

    CALL main_grid_def_init_coupling

  END SUBROUTINE main_grid_init_coupling
  ! ---------------------------------------------------------------------------
#endif
  ! ---------------------------------------------------------------------------
  SUBROUTINE main_grid_global_start

    USE messy_main_grid_def_bi, ONLY: main_grid_def_global_start

    IMPLICIT NONE

    CALL main_grid_def_global_start

  END SUBROUTINE main_grid_global_start
  ! ---------------------------------------------------------------------------

 ! ---------------------------------------------------------------------------
  SUBROUTINE main_grid_local_start

    USE messy_main_grid_def_bi, ONLY: main_grid_def_local_start

    IMPLICIT NONE

    CALL main_grid_def_local_start

  END SUBROUTINE main_grid_local_start
  ! ---------------------------------------------------------------------------
#ifdef ECHAM5
 ! ---------------------------------------------------------------------------
  SUBROUTINE main_grid_radiation

    USE messy_main_grid_def_bi, ONLY: main_grid_def_radiation

    IMPLICIT NONE

    CALL main_grid_def_radiation

  END SUBROUTINE main_grid_radiation
  ! ---------------------------------------------------------------------------
#endif
  ! ---------------------------------------------------------------------------
  SUBROUTINE main_grid_vdiff

#ifdef ICON
    USE messy_main_grid_def_bi, ONLY: main_grid_def_vdiff

    IMPLICIT NONE

    CALL main_grid_def_vdiff
#endif

  END SUBROUTINE main_grid_vdiff
  ! ---------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_free_memory

    USE messy_main_blather_bi,      ONLY: ERROR_BI
    USE messy_main_grid_def_mem_bi, ONLY: BASEGRID_ID
    USE messy_main_grid_def_bi,     ONLY: main_grid_def_free_memory
    USE messy_main_grid,            ONLY: GRID_ERROR           &
                                        , CLEAN_GEOHYBGRID_LIST
    USE messy_main_grid_trafo_scrp, ONLY: CLEAN_SCRIPDATA_LIST &
                                        , CLEAN_SCRIPGRID_LIST

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_grid_free_memory'
    INTEGER                     :: status

    CALL main_grid_def_free_memory

    CALL CLEAN_SCRIPDATA_LIST(status)
    IF (status /= 0) CALL ERROR_BI(grid_error(status), substr)
    CALL CLEAN_SCRIPGRID_LIST(status)
    IF (status /= 0) CALL ERROR_BI(grid_error(status), substr)
    CALL CLEAN_GEOHYBGRID_LIST(status)
    IF (status /= 0) CALL ERROR_BI(grid_error(status), substr)

    DEALLOCATE(BASEGRID_ID)

  END SUBROUTINE main_grid_free_memory
  !----------------------------------------------------------------------------
#ifdef ICON
  ! ---------------------------------------------------------------------------
  SUBROUTINE main_grid_set_domain(lchobj)

    USE messy_main_grid_def_bi, ONLY: main_grid_def_set_domain

    IMPLICIT NONE

    LOGICAL, INTENT(IN), OPTIONAL :: lchobj  ! reset channel object pointers

    CALL main_grid_def_set_domain(lchobj)

  END SUBROUTINE main_grid_set_domain
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------
  SUBROUTINE main_grid_set_domain_local(jb)

    USE messy_main_grid_def_bi, ONLY: main_grid_def_set_domain_local

    IMPLICIT NONE

    INTEGER, INTENT(IN)           :: jb

    CALL main_grid_def_set_domain_local(jb)

  END SUBROUTINE main_grid_set_domain_local
  ! ---------------------------------------------------------------------------
#endif
!******************************************************************************
END MODULE messy_main_grid_bi
!******************************************************************************
