MODULE mo_forcing

  USE mo_kind,       ONLY: dp, i4, i8, wp
  USE mo_parallel,   ONLY: p_pe, p_io, gather, stop_all, p_joff, global_sum
  USE mo_mpi,        ONLY: p_bcast
  USE mo_boundsexch, ONLY: bounds_exch

  USE mo_commo1,     ONLY: dt, ldtdayc, nfixYearLen,               &
                           tafo,txo,tye,fprec,fswr,                &
                           ftdew,fu10,fclou,giriv,           &
                           alat, alon, alatu, alonu, alatv, alonv, &
                           alon_g, alat_g, lbounds_exch_tp!, fslp

  USE mo_param1,     ONLY: ie, je, ie_g, je_g

  USE mo_constants,  ONLY: api, agratorad, aradtogra

  USE mo_units,      ONLY: io_stdout

  USE mo_rotation,   ONLY: rotate2_ini, rotate2_u, rotate2_v
  USE mo_grid,       ONLY: p_suchij

  USE mo_model_time, ONLY: time_desc, format_model_time, &
                           model_seconds, &
                           broadcast_time_desc, &
#ifndef __coupled
                           seconds_between_two_times, &
#endif
                           operator(+), operator(-), operator(.diff.)

  USE mo_planetary_constants, ONLY: radius, rhoref_water

#ifdef __coupled
  USE mo_couple,            ONLY: &
#  ifdef FLUXCORRECT
       couple_correct_ini, &
#  endif
       couple_get_a2o
#else
#ifdef CORE
  USE mo_ncar_ocean_fluxes, ONLY: open_core, spool_core, read_core
#else
  USE mo_omip,              ONLY: open_omip, spool_omip, read_omip, &
                                  rewind_omip
#endif
#endif


  IMPLICIT NONE

#ifndef NOCDI
  INCLUDE 'cdi.inc'
#endif

  PRIVATE

  PUBLIC :: read_namelist_forcctl, initialize_surface_forcing, &
            update_surface_forcing, finalize_surface_forcing, &
            forcing_frequency, forcing_start_time
#ifndef __coupled
  PUBLIC :: spool_forcing
#endif


  INTEGER, PARAMETER     :: maxnamelen=250  ! maximal string length
  INTEGER(i4), PARAMETER :: maxnfield=20    ! maximal number of input fields


  !----- variables that can be modified via namelist FORCCTL

  CHARACTER(LEN=20)      :: cforcdata         ! forcing data (OMIP, NCEP, ERA)
  REAL(dp)               :: forcing_frequency ! forcing intervall (sec)

  LOGICAL ::               &
#ifndef NOCDI
    lspat_interp_forcing,  & ! TRUE for bilinear HORIZONTAL interpolation of
                             ! forcing fields
    ltime_interp_forcing,  & ! TRUE for linear TIME interpolation of forcing fields
#endif
    lperiodic_forcing,     & ! TRUE to rewind forcing files at turn of year (OMIP)
                             ! (i.e. before computing january 1st)?
    ldebug_forcing,        & ! TRUE to turn on debug mode
    lwrite_forcing,        & ! TRUE to write interpolated forcing fields
    ldiff_runoff_grid        ! TRUE if runoff data are provided on other grid

#ifdef NOCDI
! op_pj_20110222+
! NOTE: the logicals are broadcasetd and therefore must not have the PARAMETER
!       attribute
!!$  LOGICAL, PARAMETER :: lspat_interp_forcing = .FALSE., &
!!$       ltime_interp_forcing = .FALSE.
  LOGICAL :: lspat_interp_forcing = .FALSE., &
       ltime_interp_forcing = .FALSE.
! op_pj_20110222-
#endif
  TYPE(time_desc), SAVE :: forcing_start_time  ! start date and time of forcing data

  !-----

  TYPE forcefile
    CHARACTER(LEN=maxnamelen) :: filename
    INTEGER  ::  &
      streamID,  & ! stream identifier
      vlistID,   & ! variable list identifier
      varID,     & ! variable identifier
      gridID,    & ! grid identifier
      taxisID,   & ! time axis identifier
      zaxisID,   & ! z-axis identifier
      current_record
  END TYPE forcefile

  TYPE field_desc                               ! filed description
    CHARACTER(LEN=maxnamelen) :: varname        ! variable name
    CHARACTER(LEN=1)   :: gtype                 ! grid type acro ('p', 'u', 'v')
    CHARACTER(LEN=20)  :: unit                  ! unit of measure
    TYPE(forcefile)    :: inp, out              ! input / output file description
    REAL(dp), POINTER  :: ft1(:,:), ft2(:,:)    ! input fields on original grid
    INTEGER(i8)        :: t1, t2                ! time
    REAL(dp), POINTER  :: dest(:,:)             ! destination
    REAL(dp), POINTER  :: spat_interp_dest(:,:) ! spatially interpolated field
    REAL(dp)           :: add, fac              ! offset and scale factor
  END TYPE field_desc

  TYPE (field_desc)    :: fdesc(maxnfield)

  TYPE grid_desc
    INTEGER(i4)           :: nlon, nlat   ! number of longitudes and latitudes
    REAL(dp), ALLOCATABLE :: lon(:,:), lat(:,:), land(:,:)
  END TYPE grid_desc

  TYPE (grid_desc), TARGET  :: fgrid   ! forcing fields grid
  TYPE (grid_desc), POINTER :: rgrid   ! runoff field grid

  TYPE ipair
    INTEGER(i4) :: a, b
  END TYPE ipair

  TYPE(ipair), ALLOCATABLE :: runoffij_cache(:,:)

  REAL(dp), ALLOCATABLE    :: area(:,:)

  INTEGER(i4) ::      &
    nread_per_day,    & ! number of forcing reads per day
    nfield=0,         & ! total number of input fields
    nfieldp=0,        & ! number of input fields to be interpolated to scalar points
    nfielduv=0,       & ! number of input fields to be interpolated to u/v-points
    indp(maxnfield),  & ! list indices of scalar fields
    induv(maxnfield), & ! list indices of vector fields
    indrunoff           ! list index of runoff field


  ! Coordinates of MPIOM grid in radiant [0,2pi], [-pi,pi]
  REAL(dp), ALLOCATABLE       :: alon1p(:,:), alat1p(:,:), &
                                 alon1u(:,:), alat1u(:,:), &
                                 alon1v(:,:), alat1v(:,:)

CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE initialize_surface_forcing

#if defined (__coupled)

#  ifdef FLUXCORRECT
    IF (p_pe == p_io) CALL couple_correct_ini
#  endif /*FLUXCORRECT*/

#else /*(__coupled)*/

   nread_per_day=NINT(86400._dp/forcing_frequency)

   IF (.NOT. lspat_interp_forcing) THEN   ! no spatial interpolation needed
#  ifdef CORE
     IF (p_pe == p_io) THEN
       WRITE(io_stdout,*) 'CALL open_core'
       CALL open_core
     ENDIF
#  else
     IF (p_pe == p_io) THEN
       WRITE(io_stdout,*) 'CALL open_omip'
       CALL open_omip
     ENDIF
#  endif /*CORE*/

    ELSE  ! spatial interpolation needed

#  ifdef NOCDI
      CALL stop_all('Error - CDI library needed to use forcing data on '// &
                     'original grid => run aborted, please recompile!')
#  else
      CALL define_forclist
      ! Get grid description of external forcing fields
      SELECT CASE (TRIM(cforcdata))
      CASE ('OMIP') ! workaround for OMIP because of different land-sea masks
        CALL get_forcing_grid(TRIM(cforcdata)//'_LSM', fgrid, 1.e-4_dp)
        ALLOCATE(rgrid)
        CALL get_forcing_grid(TRIM(cforcdata)//'_LSM', rgrid)
      CASE DEFAULT
        CALL get_forcing_grid(TRIM(cforcdata)//'_LSM', fgrid)
        IF (ldiff_runoff_grid) THEN
          ALLOCATE(rgrid)
          CALL get_forcing_grid(TRIM(cforcdata)//'_LSM_RUNOFF', rgrid)
        ELSE
          rgrid => fgrid
        ENDIF
      END SELECT

      CALL model_grid_coordinates

      CALL rotate2_ini


      CALL runoff_ini(rgrid%nlon, rgrid%nlat, rgrid%lat, rgrid%lon, rgrid%land)

      IF (lwrite_forcing) CALL write_forcing_ini

#  endif /*NOCDI*/

    ENDIF

#endif /*(__coupled)*/

  END SUBROUTINE initialize_surface_forcing

!------------------------------------------------------------------------------

  SUBROUTINE update_surface_forcing(model_time, model_start_time, ldtrun)
    TYPE(time_desc), INTENT(IN) :: model_time, model_start_time
    INTEGER, INTENT(IN) :: ldtrun
#ifndef __coupled
    TYPE(time_desc) :: forc_time, forc_time_out
    REAL(dp) :: wgt1, wgt2 , forc_period_fraction
    REAL(dp), POINTER :: temp(:,:)
    INTEGER :: i, vdate, vtime, req_idx
    INTEGER(i8) :: tc, mod_sec, forc_sec
#endif

#if defined (__coupled)
    ! Get data from coupler
    CALL couple_get_a2o(ldtrun)
#else /*(__coupled)*/


    IF (.NOT. lspat_interp_forcing) THEN   ! no spatial interpolation needed
      IF (MOD(ldtdayc-1,NINT(forcing_frequency/dt)) .EQ. 0) THEN
#  ifdef CORE
        WRITE(io_stdout,*) 'CALL read_core'
        CALL read_core
#  else
        WRITE(io_stdout,*) 'CALL read_omip'
        CALL read_omip
#  endif
      ENDIF

    ELSE   ! spatial interpolation needed

      forc_time = forcing_start_time + seconds_between_two_times(model_time,model_start_time)
      IF (lwrite_forcing) forc_time_out = forc_time

      IF (ltime_interp_forcing) THEN
        tc = next_forcing_time_code(forc_time, forcing_start_time%year)
        forc_sec = tc * INT(forcing_frequency,i8)

        !  coefficients for linear time interpolation
        forc_period_fraction = (REAL(model_seconds(forc_time,forcing_start_time%year),dp) &
                              - REAL(forc_sec,dp)) / forcing_frequency

        wgt2=forc_period_fraction + REAL(FLOOR(1._dp-forc_period_fraction),dp)
        wgt1=1._dp-wgt2

        mod_sec = model_seconds(model_time,model_start_time%year)
        forc_time = forc_time + (forc_sec - mod_sec)
      ENDIF
      req_idx=record_idx(forc_time)
      DO i=1, nfield
        ! Inquire the time step
        CALL prepare_record_in(fdesc(i), forc_time, req_idx)
      ENDDO
      IF (ltime_interp_forcing) THEN
        IF (tc > fdesc(1)%t2) THEN
          DO i=1, nfield
            temp => fdesc(i)%ft1
            fdesc(i)%ft1 => fdesc(i)%ft2
            fdesc(i)%t1=fdesc(i)%t2
            fdesc(i)%ft2 => temp
            fdesc(i)%t2=tc
            fdesc(i)%spat_interp_dest => temp
          ENDDO
          CALL read_interpolate_forcing(forc_time, req_idx)
        ENDIF

        IF (ldebug_forcing .AND. p_pe == p_io) &
            WRITE (0,'(a,3i8,2f7.4)') 'Time codes and weighting factors', &
                            fdesc(1)%t1, fdesc(1)%t2, tc, wgt1, wgt2

        DO i=1,nfield
          CALL time_interp_forcing(wgt1, wgt2, fdesc(i))
        ENDDO
      ELSE
        IF (MOD(ldtdayc-1,NINT(forcing_frequency/dt)) .EQ. 0) THEN
          CALL read_interpolate_forcing(forc_time, req_idx)
        ENDIF
      ENDIF

      IF (lwrite_forcing) THEN
        vdate=forc_time_out%year*10000 + forc_time_out%month*100 + forc_time_out%mday
        vtime=forc_time_out%hour*10000 + forc_time_out%minute*100 + forc_time_out%second
        IF (ltime_interp_forcing .OR. vtime == 0) THEN
          DO i=1,nfield
            CALL write_forcing(fdesc(i), vdate, vtime)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

#endif /*(__coupled)*/
  END SUBROUTINE update_surface_forcing

!------------------------------------------------------------------------------

  SUBROUTINE finalize_surface_forcing
    INTEGER(i4) :: ii

      DO ii=1,nfield
        IF (ltime_interp_forcing) DEALLOCATE(fdesc(ii)%ft1, fdesc(ii)%ft2)
        IF (p_pe == p_io) THEN
          CALL close_forcing_file(fdesc(ii)%inp%filename,fdesc(ii)%inp%streamID)
          IF (lwrite_forcing) CALL close_forcing_file(fdesc(ii)%out%filename,fdesc(ii)%out%streamID)
          IF (lwrite_forcing) CALL destroy_objects(fdesc(ii))
        ENDIF
      ENDDO
  END SUBROUTINE finalize_surface_forcing

!------------------------------------------------------------------------------

  !> Define and read namelist FORCCTL
  SUBROUTINE read_namelist_forcctl(model_start_time, io_in_forcctl, ierror)

    TYPE(time_desc), INTENT(IN) :: model_start_time
    INTEGER, INTENT(in) :: io_in_forcctl
    INTEGER, INTENT(OUT) :: ierror

    NAMELIST /forcctl/ cforcdata, forcing_frequency, lwrite_forcing, &
#ifndef NOCDI
                       ltime_interp_forcing, lspat_interp_forcing, &
#endif
                       lperiodic_forcing, &
                       ldebug_forcing, ldiff_runoff_grid, forcing_start_time

    ! Set default values
    cforcdata = 'OMIP'              ! forcing data set
    forcing_frequency = 86400._dp   ! forcing data period (sec)
    lwrite_forcing = .FALSE.        ! write interpolated forcing fields
#ifndef NOCDI
    ltime_interp_forcing = .FALSE.  ! perform time interpolation
    lspat_interp_forcing = .FALSE.  ! perform spatial interpolation
#endif
    lperiodic_forcing = .FALSE.     ! reposition forcing files to initial point at turn of year
    ldebug_forcing = .FALSE.        ! turn on debug mode
    ldiff_runoff_grid = .TRUE.      ! differing original grid for runoff data
    forcing_start_time = model_start_time     ! forcing start year

    ! Read namelist forcctl
    IF ( p_pe == p_io ) READ(io_in_forcctl,forcctl,iostat=ierror)

    ! Broadcast external namelist settings
    CALL p_bcast(cforcdata, p_io)
    CALL p_bcast(forcing_frequency, p_io)
    CALL p_bcast(lwrite_forcing, p_io)
    CALL p_bcast(ltime_interp_forcing, p_io)
    CALL p_bcast(lspat_interp_forcing, p_io)
    CALL p_bcast(lperiodic_forcing, p_io)
    CALL p_bcast(ldebug_forcing, p_io)
    CALL p_bcast(ldiff_runoff_grid, p_io)
    CALL broadcast_time_desc(forcing_start_time)

    SELECT CASE(TRIM(cforcdata))
      CASE('NCEP','ERA','OMIP')

      CASE DEFAULT
        CALL stop_all('Error - wrong forcing data type: '//TRIM(cforcdata)// &
                      ' => run aborted')
    END SELECT

    IF ( MOD(forcing_frequency,dt) .GT. 0._dp ) THEN
      CALL stop_all('STOP : timestep does not match with forcing frequency ')
    ENDIF

#ifndef NOCDI
    IF (.NOT. lspat_interp_forcing) ltime_interp_forcing=.FALSE.
#endif
    IF (.NOT. lspat_interp_forcing) lwrite_forcing=.FALSE.

  END SUBROUTINE read_namelist_forcctl

!------------------------------------------------------------------------------

  SUBROUTINE define_forclist

    ! Purpose: define forcing fields
    ! (variable name, target grid acro, unit of measure, destination field)

    CALL forclist_add('TEM'   ,'p'  ,'K'      , tafo)
    CALL forclist_add('TDEW'  ,'p'  ,'K'      , ftdew)
    CALL forclist_add('PREC'  ,'p'  ,'m s-1'  , fprec)
    CALL forclist_add('SWRAD' ,'p'  ,'W m-2'  , fswr)
    CALL forclist_add('WIND10','p'  ,'m s-1'  , fu10)
    CALL forclist_add('CLOUD' ,'p'  ,' '      , fclou)
    CALL forclist_add('WIX'   ,'u'  ,'Pa'     , txo)
    CALL forclist_add('WIY'   ,'v'  ,'Pa'     , tye)
   !CALL forclist_add('PRESS' ,'p'  ,'hPa'    , fslp)
   !CALL forclist_add('LWRAD' ,'p'  ,'W m-2'  , xxx)
   !CALL forclist_add('U10'   ,'u'  ,'m s-1'  , xxx)
   !CALL forclist_add('V10'   ,'v'  ,'m s-1'  , xxx)
    CALL forclist_add('RIV'   ,'p'  ,'m s-1'  , giriv)

  END SUBROUTINE define_forclist

!------------------------------------------------------------------------------

  !> Fill forcing field description structure
  SUBROUTINE forclist_add(varname,gtype,unit,dest)

    CHARACTER(LEN=*),  INTENT(IN)   :: varname
    CHARACTER(LEN=*),  INTENT(IN)   :: gtype
    CHARACTER(LEN=*),  INTENT(IN)   :: unit
    REAL(dp), INTENT(INOUT), TARGET :: dest(:,:)

    nfield=nfield+1
    IF (TRIM(varname) .EQ. 'RIV') THEN
      indrunoff=nfield
    ELSEIF (gtype .EQ. 'p') THEN
      nfieldp=nfieldp+1
      indp(nfieldp)=nfield
    ELSEIF (gtype .EQ. 'u' .OR. gtype .EQ. 'v') THEN
      nfielduv=nfielduv+1
      induv(nfielduv)=nfield
    ENDIF

    fdesc(nfield)%varname=varname
    fdesc(nfield)%gtype=gtype
    fdesc(nfield)%unit=unit
    fdesc(nfield)%inp%current_record = 0
    fdesc(nfield)%out%current_record = 0
    fdesc(nfield)%inp%filename = ''
    fdesc(nfield)%out%filename = ''
    fdesc(nfield)%t1 = -1_i8
    fdesc(nfield)%t2 = -1_i8
    fdesc(nfield)%dest => dest
    fdesc(nfield)%spat_interp_dest => dest
    if (ltime_interp_forcing) then
      ALLOCATE(fdesc(nfield)%ft1(ie, je), fdesc(nfield)%ft2(ie, je))
    end if

    ! Define data dependent offset and scaling factor
    fdesc(nfield)%add = 0.0_dp
    fdesc(nfield)%fac = 1.0_dp
    SELECT CASE (TRIM(fdesc(nfield)%varname))
      CASE ('TEM')
        fdesc(nfield)%add = -273.15_dp    ! convert to degC
      CASE ('CLOUD')
        IF (TRIM(cforcdata) .EQ. 'NCEP') fdesc(nfield)%fac = 1.0_dp/100.0_dp
      CASE ('PREC')
        IF (TRIM(cforcdata) .EQ. 'NCEP') fdesc(nfield)%fac = 1.0_dp/1000.0_dp
      CASE ('U10', 'V10')
        IF (TRIM(cforcdata) .EQ. 'NCEP') fdesc(nfield)%fac = -1.0_dp
      CASE ('WIX', 'WIY')
        IF (TRIM(cforcdata) .EQ. 'NCEP') fdesc(nfield)%fac = -1.0_dp/rhoref_water
        IF (TRIM(cforcdata) .EQ. 'OMIP') fdesc(nfield)%fac =  1.0_dp/rhoref_water
        IF (TRIM(cforcdata) .EQ. 'ERA' ) fdesc(nfield)%fac =  1.0_dp/rhoref_water
    END SELECT

  END SUBROUTINE forclist_add

!------------------------------------------------------------------------------

  SUBROUTINE get_forcing_grid(filename, grid, land_frac)

    CHARACTER(LEN=*), INTENT(IN)    :: filename   ! filename for land sea mask
                                                  ! (self-described file format)
    TYPE(grid_desc), INTENT(OUT)    :: grid

    REAL(dp), INTENT(IN), OPTIONAL  :: land_frac

#ifndef NOCDI
    REAL(dp), ALLOCATABLE           :: lon(:), lat(:), lsm(:,:)

    INTEGER(i4)                     :: streamID, vlistID, &
                                       gridID, nlon, nlat, nmiss, i, j
    IF (p_pe==p_io) THEN

      ! Get land-sea mask and grid coordinates

      WRITE(0,*) 'Open input file '//TRIM(filename)
      streamID = streamOpenRead(TRIM(filename))
      IF ( streamID < 0 ) THEN
        CALL stop_all ('Problem opening file '//TRIM(filename)//': '  &
                       //cdiStringError(streamID))
      ENDIF

      ! Get variable list identifier
      vlistID = streamInqVlist(streamID)

      ! Get grid identifier
      gridID = vlistInqVarGrid(vlistID, 0)

      ! Get grid dimensions nlon and nlat
      nlon=gridInqXSize(gridID)
      nlat=gridInqYSize(gridID)

    ENDIF

    CALL p_bcast(nlon,p_io)
    CALL p_bcast(nlat,p_io)

    ALLOCATE(lon(nlon), lat(nlat), lsm(nlon,nlat))
    ALLOCATE(grid%lon(nlon+2,nlat+2), grid%lat(nlon+2,nlat+2), grid%land(nlon+2,nlat+2))

    IF (p_pe==p_io) THEN

      ! Get grid coordinates lon and lat
      nlon=gridInqXvals(gridID,lon)
      nlat=gridInqYvals(gridID,lat)
      IF(ldebug_forcing) WRITE(0,'(a20,g14.6,a,g14.6)') 'FORCING FIELD LONs:',lon(1),'...',lon(nlon)
      IF(ldebug_forcing) WRITE(0,'(a20,g14.6,a,g14.6)') 'FORCING FIELD LATs:',lat(1),'...',lat(nlat)

      ! Read land sea mask
      CALL streamReadVar(streamID, 0, lsm(:,:), nmiss )
      IF (ldebug_forcing) WRITE(0,'(a20,3g11.3)') 'FORCING FIELD LSM', minval(lsm), &
                                      maxval(lsm), sum(lsm)/REAL(nlon*nlat,dp)
    ENDIF

    CALL p_bcast(lon,p_io)
    CALL p_bcast(lat,p_io)
    CALL p_bcast(lsm,p_io)

    ! Construct 2d longutide and latitude arrays
    grid%nlon=nlon
    grid%nlat=nlat
    DO i=1,nlon
      DO j=1,nlat
        grid%lon(i+1,j+1) = lon(i)
        grid%lat(i+1,j+1) = lat(j)
      ENDDO
    ENDDO

    ! cyclic east-west boundary
    DO j=2,nlat+1
      grid%lon(1,j)=grid%lon(nlon+1,j) - 360._dp
      grid%lat(1,j)=grid%lat(nlon+1,j)
      grid%lon(nlon+2,j)=grid%lon(2,j) + 360._dp
      grid%lat(nlon+2,j)=grid%lat(2,j)
    ENDDO

    ! "extrapolate" coordinates for North Pole and South Pole
    grid%lon(:,1)=grid%lon(:,2)
    grid%lat(:,1)=180._dp - grid%lat(:,2)
    grid%lon(:,nlat+2)=grid%lon(:,nlat+1)
    grid%lat(:,nlat+2)=-180._dp - grid%lat(:,nlat+1)

    ! convert to radiant
    grid%lon=grid%lon*agratorad
    grid%lat=grid%lat*agratorad

    ! Construct land array

    IF (PRESENT(land_frac)) THEN
      grid%land(2:nlon+1,2:nlat+1) = MERGE(1._dp, lsm, lsm .GT. land_frac)
    ELSE
      grid%land(2:nlon+1,2:nlat+1) = lsm
    ENDIF
    grid%land(:,1)=0._dp
    grid%land(:,nlat+2)=1._dp
    grid%land(1,:)=grid%land(nlon+1,:)
    grid%land(nlon+2,:)=grid%land(2,:)

    IF (p_pe == p_io) CALL streamClose(streamID)

#endif
  END SUBROUTINE get_forcing_grid

!------------------------------------------------------------------------------

  SUBROUTINE model_grid_coordinates

    ! Convert MPIOM grid coordinates to radiant [0,2pi],[-pi,pi]

    ALLOCATE(alon1p(ie,je), alat1p(ie,je), alon1u(ie,je), alat1u(ie,je), &
             alon1v(ie,je), alat1v(ie,je))

    ! scalar points
    alat1p=alat*agratorad
    alon1p=alon*agratorad
    WHERE (alon1p .GE. (2._dp*api))
      alon1p=alon1p-2._dp*api
    ENDWHERE
    WHERE (alon1p .LT. 0._dp)
      alon1p=alon1p+2._dp*api
    ENDWHERE

    ! u-points
    alat1u=alatu*agratorad
    alon1u=alonu*agratorad
    WHERE (alon1u .GE. (2._dp*api))
      alon1u=alon1u-2._dp*api
    ENDWHERE
    WHERE (alon1u .LT. 0._dp)
      alon1u=alon1u+2._dp*api
    ENDWHERE

    ! v-points
    alat1v=alatv*agratorad
    alon1v=alonv*agratorad
    WHERE (alon1v .GE. (2._dp*api))
      alon1v=alon1v-2._dp*api
    ENDWHERE
    WHERE (alon1v .LT. 0._dp)
      alon1v=alon1v+2._dp*api
    ENDWHERE

  END SUBROUTINE model_grid_coordinates

!------------------------------------------------------------------------------

  SUBROUTINE open_forcing_file(fdesc, newfn)

    TYPE(field_desc),          INTENT(INOUT) :: fdesc
    CHARACTER(LEN=maxnamelen), INTENT(IN) :: newfn

#ifndef NOCDI

    ! Open input file
    WRITE(0,*) 'Open input  file '//TRIM(newfn)
    fdesc%inp%streamID = streamOpenRead(TRIM(newfn))
    IF ( fdesc%inp%streamID < 0 ) THEN
       CALL stop_all ('Problem opening file '//TRIM(newfn)//': '  &
            //cdiStringError(fdesc%inp%streamID))
    ENDIF

    ! Get the variable list identifier
    fdesc%inp%vlistID = streamInqVlist(fdesc%inp%streamID)

    ! Set the variable IDs
    fdesc%inp%varID = 0

    ! Get the grid ID
    fdesc%inp%gridID = vlistInqVarGrid(fdesc%inp%vlistID, fdesc%inp%varID)

    ! Get the time axis from the variable list
    fdesc%inp%taxisID = vlistInqTaxis(fdesc%inp%vlistID)

    ! after opening file pointer is at first record
    fdesc%inp%current_record = 0

#endif
  END SUBROUTINE open_forcing_file

!------------------------------------------------------------------------------

  SUBROUTINE close_forcing_file(filename,streamID)
    CHARACTER(LEN=*), INTENT(INOUT) :: filename
    INTEGER, INTENT(IN) :: streamID
#ifndef NOCDI
    CALL streamClose(streamID)
    filename = ''
#endif
  END SUBROUTINE close_forcing_file

!------------------------------------------------------------------------------

  SUBROUTINE prepare_record_in(fdesc, req_time, req_idx)

    TYPE(field_desc), INTENT(INOUT) :: fdesc
    TYPE(time_desc),  INTENT(IN)    :: req_time
    INTEGER,          INTENT(IN)    :: req_idx
#ifndef NOCDI
    INTEGER(i4)   :: spool_idx, nrec
    CHARACTER(LEN=maxnamelen) :: new_fname

    IF (p_pe == p_io) THEN
      ! open file if necessary
      new_fname = spec_filename(fdesc%varname, req_time%year)
      IF (fdesc%inp%filename .ne. new_fname .OR. &      ! turn of year reached
        (TRIM(cforcdata) .EQ. 'OMIP' .AND. req_idx .EQ. 0)) THEN
        IF (TRIM(fdesc%inp%filename) .ne. '') THEN
          CALL close_forcing_file(fdesc%inp%filename,fdesc%inp%streamID)
        END IF
        fdesc%inp%filename = new_fname
        CALL open_forcing_file(fdesc, new_fname)
      END IF
      ! Spool to the requested time step
      spool_idx = fdesc%inp%current_record
      IF (spool_idx > req_idx) THEN
        CALL stop_all('Unhandled')
      END IF
      DO WHILE (spool_idx < req_idx)
        nrec = streamInqTimestep(fdesc%inp%streamID, spool_idx)
        spool_idx=spool_idx+1
        IF (nrec == 0 ) CALL stop_all('Error in prepare_record_in')
      ENDDO
      fdesc%inp%current_record = spool_idx
    ENDIF

#endif
  END SUBROUTINE prepare_record_in

!------------------------------------------------------------------------------

#ifndef __coupled

  SUBROUTINE spool_forcing

    IF (.NOT. lspat_interp_forcing) THEN
#ifdef CORE
      WRITE(io_stdout,*) 'CALL spool_core'
      IF (lperiodic_forcing) &
        CALL stop_all('periodic forcing not implemented for core')
      CALL spool_core(nread_per_day)
#else
      WRITE(io_stdout,*) 'CALL spool_omip'
      IF (lperiodic_forcing .AND. p_pe == p_io) CALL rewind_omip
      CALL spool_omip(nread_per_day)
#endif
     ENDIF

  END SUBROUTINE spool_forcing

#endif/*ndef __coupled */

!------------------------------------------------------------------------------

  SUBROUTINE read_interpolate_forcing(req_time, req_idx)
    TYPE(time_desc), INTENT(IN) :: req_time
    INTEGER,         INTENT(IN) :: req_idx
    INTEGER(i4)              :: ii, ip, iuv
    REAL(dp), ALLOCATABLE    :: pfld(:,:,:), prunoff(:,:), uvfld(:,:,:),  &
                                pfld_interp(:,:,:), prunoff_interp(:,:), uvfld_interp(:,:,:)


    ALLOCATE(pfld(fgrid%nlon+2,fgrid%nlat+2,nfieldp), &
             pfld_interp(ie,je,nfieldp), &
             prunoff(rgrid%nlon,rgrid%nlat), &
             prunoff_interp(ie,je), &
             uvfld(fgrid%nlon+2,fgrid%nlat+2,nfielduv),  &
             uvfld_interp(ie,je,nfielduv))

    IF (p_pe == p_io) &
      WRITE(0,'(a,i8)') TRIM(cforcdata)//': read forcing:  '// &
                        TRIM(format_model_time(req_time))//  &
                        '  record from file:', req_idx+1
    ip=0
    iuv=0
    DO ii=1,nfield
      IF (TRIM(fdesc(ii)%varname) .EQ. "RIV") THEN
        CALL read_one_forcing_stream(fdesc(ii), req_idx, &
                                    rgrid%nlon, rgrid%nlat, prunoff, .FALSE.)
      ELSEIF (TRIM(fdesc(ii)%gtype) .EQ. 'p') THEN
        ip=ip+1
        CALL read_one_forcing_stream(fdesc(ii), req_idx, &
                                    fgrid%nlon+2, fgrid%nlat+2, pfld(:,:,ip))

      ELSEIF (TRIM(fdesc(ii)%gtype) .EQ. 'u'                      &
         .OR. TRIM(fdesc(ii)%gtype) .EQ. 'v') THEN
        iuv=iuv+1
        CALL read_one_forcing_stream(fdesc(ii), req_idx, &
                                    fgrid%nlon+2, fgrid%nlat+2, uvfld(:,:,iuv))
      ENDIF
    ENDDO

    CALL p_bcast(prunoff, p_io)
    CALL p_bcast(pfld, p_io)
    CALL p_bcast(uvfld,p_io)

    ! Interpolation to scalar p-grid

    CALL interpolate_bln(nfieldp,fgrid%nlon,fgrid%nlat,ie,je,fgrid%lat,fgrid%lon,alat1p,alon1p, &
                                 pfld,pfld_interp,fgrid%land)

    CALL runoff(rgrid%nlon, rgrid%nlat, ie, je, prunoff, prunoff_interp, &
                rgrid%land(2:rgrid%nlon+1,2:rgrid%nlat+1))

    CALL bounds_exch(1,'p',pfld_interp,'mo_forcing p')
    CALL bounds_exch(1,'p',prunoff_interp,'mo_forcing p')

    DO ii=1,nfieldp
      ip=indp(ii)
      fdesc(ip)%spat_interp_dest=pfld_interp(:,:,ii)
    ENDDO

    fdesc(indrunoff)%spat_interp_dest=prunoff_interp

    ! Interpolate to u-grid
    CALL interpolate_bln(nfielduv,fgrid%nlon,fgrid%nlat,ie,je,fgrid%lat,fgrid%lon,alat1u,alon1u, &
                                 uvfld,uvfld_interp,fgrid%land)

    ! rotate vectors
    DO ii=1,nfielduv,2
      iuv=induv(ii)
      CALL rotate2_u(uvfld_interp(:,:,ii),uvfld_interp(:,:,ii+1),ie,je)
      CALL bounds_exch(1,'u',uvfld_interp(:,:,ii),'mo_forcing u')

      fdesc(iuv)%spat_interp_dest=uvfld_interp(:,:,ii)
    ENDDO

    ! Interpolate to v-grid
    CALL interpolate_bln(nfielduv,fgrid%nlon,fgrid%nlat,ie,je,fgrid%lat,fgrid%lon,alat1v,alon1v, &
                                 uvfld,uvfld_interp,fgrid%land)

    ! rotate vectors
    DO ii=1,nfielduv,2
      iuv=induv(ii+1)
      CALL rotate2_v(uvfld_interp(:,:,ii),uvfld_interp(:,:,ii+1),ie,je)
      CALL bounds_exch(1,'v',uvfld_interp(:,:,ii+1),'mo_forcing v')

      fdesc(iuv)%spat_interp_dest=uvfld_interp(:,:,ii+1)
    ENDDO

  END SUBROUTINE read_interpolate_forcing

!------------------------------------------------------------------------------

  SUBROUTINE read_one_forcing_stream(fdesc, req_idx, nx, ny, ffield, lglobhalo)

    TYPE(field_desc), INTENT(INOUT) :: fdesc
    INTEGER,          INTENT(IN) :: req_idx
    INTEGER(i4),      INTENT(IN) :: nx, ny
    LOGICAL,          INTENT(IN), OPTIONAL :: lglobhalo

    REAL(dp),         INTENT(OUT) :: ffield(nx,ny)
#ifndef NOCDI
    INTEGER  :: vdate, vtime
    INTEGER  :: nmiss, nrec
    LOGICAL  :: lglh

    IF (p_pe == p_io) THEN
      nrec = streamInqTimestep(fdesc%inp%streamID, req_idx)
      IF (nrec == 0 ) &
        CALL stop_all('Error - requested time step not available in '// &
                      'forcing file '//fdesc%inp%filename)

      ! Get verification date and time
      vdate = taxisInqVdate(fdesc%inp%taxisID)
      vtime = taxisInqVtime(fdesc%inp%taxisID)
      IF (ldebug_forcing) WRITE(0,'(a,i8.8,a,i6.6)') 'VDATE: ', vdate, ' VTIME: ', vtime

      ! Read input field
      lglh=.TRUE.
      IF (PRESENT(lglobhalo)) lglh=lglobhalo

      IF (lglh) THEN
        CALL streamReadVar(fdesc%inp%streamID, fdesc%inp%varID, &
             ffield(2:nx-1,2:ny-1), nmiss)

        ! cyclic east-west boundary
        ffield(1,:)=ffield(nx-1,:)
        ffield(nx,:)=ffield(2,:)

        ! treatment of North and South Poles
        ffield(:,1)=ffield(:,2)
        ffield(:,ny)=ffield(:,ny-1)

        IF (ldebug_forcing) WRITE(0,'(a,3g14.6)') TRIM(fdesc%varname)//': ', &
                            minval(ffield(2:nx-1,2:ny-1)), maxval(ffield(2:nx-1,2:ny-1)), &
                            sum(ffield(2:nx-1,2:ny-1))/REAL((nx-2)*(ny-2),dp)
      ELSE
        CALL streamReadVar(fdesc%inp%streamID, fdesc%inp%varID, ffield, nmiss)
        IF (ldebug_forcing) WRITE(0,'(a,3g14.6)') TRIM(fdesc%varname)//': ', &
                            minval(ffield), maxval(ffield), sum(ffield)/REAL(nx*ny,dp)
      ENDIF


      ! Modify input data if necessary
      SELECT CASE (TRIM(fdesc%varname))
        CASE ('PREC', 'SWRAD', 'LWRAD')
          ffield=MERGE(0._dp, ffield, ffield .LT. 0._dp)
      END SELECT
      ffield = ffield * fdesc%fac + fdesc%add

    ENDIF
#endif
  END SUBROUTINE read_one_forcing_stream

!------------------------------------------------------------------------------

  SUBROUTINE interpolate_bln(numfi,nlon,nlat,ie,je,alat2,alon2,alat1,alon1,&
                             soe,sreg,land)

  ! bilinear interpolation from one grid(nlon,nlat) to another grid(ie,je)

    INTEGER(i4), INTENT(IN)     :: numfi, nlon, nlat
    INTEGER(i4), INTENT(IN)     :: ie, je
    REAL(dp),    INTENT(IN)     :: alat2(nlon+2, nlat+2), alon2(nlon+2, nlat+2)
    REAL(dp),    INTENT(IN)     :: alat1(ie, je), alon1(ie, je)
    REAL(dp),    INTENT(INOUT)  :: soe(nlon+2, nlat+2, numfi)
    REAL(dp),    INTENT(OUT)    :: sreg(ie, je, numfi)
    REAL(dp),    INTENT(IN)     :: land(nlon+2, nlat+2)

    INTEGER(i4)                 :: i, j, l, m, n, iter, il, ir, jo, ju, jlu, ilu, nland
    REAL(dp), PARAMETER         :: epsilon=1.e-12_dp
    REAL(dp)                    :: g(nlon+2,nlat+2), hg(nlon+2,nlat+2), &
                                   hx(nlon+2,nlat+2,numfi), hhx(nlon+2,nlat+2,numfi), &
                                   rsumg, alpha, beta, wwwalpha, wwwbeta

    ! diffusion into land
    hhx=soe

    WHERE (land .GE. 0.5_dp)
      hg=epsilon
    ELSEWHERE
      hg=1._dp
    ENDWHERE

    DO iter=1,100
      g=hg
      hx=hhx
      DO j=1,nlat+2
        jo=max(j-1,1)
        ju=min(j+1,nlat+2)
        DO i=1,nlon+2
          il=i-1
          IF(il.lt.1) il=nlon
          ir=i+1
          IF(ir.gt.nlon+2) ir=3
          rsumg = 0._wp
          IF (land(i,j) .GE. 0.5_dp) THEN
            rsumg = (4._wp * g(i, j) + g(il, j) + g(ir, j) &
                 + g(i, jo) + g(i, ju)) / 8._wp

            hg(i,j) = MIN(rsumg,0.125_wp)

            DO l=1,numfi
              hhx(i,j,l) = (4._wp * hx(i, j, l) * g(i, j) &
                   + hx(il, j, l) * g(il, j) &
                   + hx(ir, j, l) * g(ir, j) &
                   + hx(i, jo, l) * g(i, jo) &
                   + hx(i, ju, l) * g(i, ju)) / 8._wp

              hhx(i,j,l)=hhx(i,j,l)/rsumg
            ENDDO
          ENDIF
        ENDDO
      ENDDO

      nland=0
      DO i=1,nlon
        DO j=1,nlat
          IF (hg(i,j) .LE. 2._wp * epsilon) nland=nland+1
        ENDDO
      ENDDO
      !WRITE(0,*) iter,' nland: ', nland
    ENDDO

    DO j=1,nlat+2
      DO i=1,nlon+2
        g(i,j)=hg(i,j)
        DO l=1,numfi
          soe(i,j,l)=hhx(i,j,l)
          IF (ABS(soe(i,j,l)) .LT. epsilon) soe(i, j, l) = 0._wp
        ENDDO
      ENDDO
    ENDDO

    DO m=1,ie
      DO n=1,je

      !        point left above
        DO j=1,nlat+2
          IF (alat2(2, j) .GE. alat1(m, n)) jo=j
        ENDDO
        DO i=1,nlon+2
          IF (alon2(i, 2) .LE. alon1(m, n)) il=i
        ENDDO

      !       point left above --> left below

        jlu=jo+1
        ilu=il

        wwwalpha=alat2(ilu,jlu)-alat1(m,n)
        wwwbeta=alon2(ilu,jlu)-alon1(m,n)
        IF (wwwbeta .GE. api) wwwbeta = wwwbeta - 2._wp * api
        IF (wwwbeta .LE. -api) wwwbeta = wwwbeta + 2._wp * api

        alpha=(wwwalpha)/(alat2(ilu,jlu)-alat2(ilu,jlu-1))
        beta=(wwwbeta)/(alon2(ilu,jlu)-alon2(ilu+1,jlu))
        DO i=1,numfi
          sreg(m,n,i)=alpha*beta*soe(ilu+1,jlu-1,i)*g(ilu+1,jlu-1)   &
               + (1._wp - alpha) * (1._wp - beta) * soe(ilu, jlu, i) * g(ilu, jlu) &
               + (1._wp - alpha) * beta * soe(ilu+1, jlu, i) * g(ilu+1, jlu) &
               + alpha * (1._wp - beta) * soe(ilu, jlu-1, i) * g(ilu, jlu-1)
          sreg(m,n,i)=sreg(m,n,i)/ &
               (alpha*beta*g(ilu+1,jlu-1) &
               + (1._wp - alpha) * (1._wp - beta) * g(ilu, jlu) &
               + (1._wp - alpha) * beta * g(ilu+1, jlu) &
               + alpha * (1._wp - beta) * g(ilu, jlu-1))
        ENDDO

      ENDDO
    ENDDO

    RETURN

  END SUBROUTINE interpolate_bln

!------------------------------------------------------------------------------

  SUBROUTINE write_forcing_ini
#ifndef NOCDI

    INTEGER(i4)           :: p_gridID, u_gridID, v_gridID, zaxisID, ii
    REAL(dp), ALLOCATABLE :: alonu_g(:,:), alatu_g(:,:), &
                             alonv_g(:,:), alatv_g(:,:)
    INTEGER               :: rdate, rtime

    IF ( p_pe == p_io ) THEN
      ALLOCATE(alonu_g(ie_g,je_g))
      ALLOCATE(alatu_g(ie_g,je_g))
      ALLOCATE(alonv_g(ie_g,je_g))
      ALLOCATE(alatv_g(ie_g,je_g))
    ELSE
      ALLOCATE(alonu_g(0,0))
      ALLOCATE(alatu_g(0,0))
      ALLOCATE(alonv_g(0,0))
      ALLOCATE(alatv_g(0,0))
    ENDIF

    CALL gather(alonu,alonu_g,p_io)
    CALL gather(alatu,alatu_g,p_io)
    CALL gather(alonv,alonv_g,p_io)
    CALL gather(alatv,alatv_g,p_io)


    ! Create grids
    IF ( p_pe == p_io ) THEN
      ! p-grid
      p_gridID = gridCreate(grid_curvilinear, (ie_g*je_g))
      CALL gridDefXsize(p_gridID, ie_g)
      CALL gridDefYsize(p_gridID, je_g)
      CALL gridDefXvals(p_gridID, alon_g)
      CALL gridDefYvals(p_gridID, alat_g)

      ! u-grid
      u_gridID = gridCreate(grid_curvilinear, (ie_g*je_g))
      CALL gridDefXsize(u_gridID, ie_g)
      CALL gridDefYsize(u_gridID, je_g)
      CALL gridDefXvals(u_gridID, alonu_g)
      CALL gridDefYvals(u_gridID, alatu_g)

      ! v-grid
      v_gridID = gridCreate(grid_curvilinear, (ie_g*je_g))
      CALL gridDefXsize(v_gridID, ie_g)
      CALL gridDefYsize(v_gridID, je_g)
      CALL gridDefXvals(v_gridID, alonv_g)
      CALL gridDefYvals(v_gridID, alatv_g)


      ! Create z-axis
      zaxisID = zaxisCreate(zaxis_surface, 1)
      CALL zaxisDefLevels(zaxisID, (/0.0_dp/))

      DO ii=1,nfield
        ! Create variable list
        fdesc(ii)%out%vlistID = vlistCreate()

        ! Define grids and variables
        fdesc(ii)%out%zaxisID = zaxisID
        SELECT CASE (TRIM(fdesc(ii)%gtype))
        CASE ('p')
          fdesc(ii)%out%gridID = p_gridID
        CASE ('u')
          fdesc(ii)%out%gridID = u_gridID
        CASE ('v')
          fdesc(ii)%out%gridID = v_gridID
        END SELECT

        fdesc(ii)%out%varID = vlistDefVar(fdesc(ii)%out%vlistID, &
             fdesc(ii)%out%gridID, fdesc(ii)%out%zaxisID, time_variable)

        ! Define variable names
        CALL vlistDefVarName(fdesc(ii)%out%vlistID, fdesc(ii)%out%varID, &
             fdesc(ii)%varname)

        ! Create a time axis, define calendar and reference date/time
        fdesc(ii)%out%taxisID = taxisCreate(taxis_relative)
        IF (nfixYearLen .EQ. 360) THEN
          CALL taxisDefCalendar(fdesc(ii)%out%taxisID, calendar_360days)
        ELSEIF (nfixYearLen .EQ. 365) THEN
          CALL taxisDefCalendar(fdesc(ii)%out%taxisID, calendar_365days)
        ELSE
          CALL taxisDefCalendar(fdesc(ii)%out%taxisID, calendar_proleptic)
        ENDIF

        rdate = forcing_start_time%year * 10000 + forcing_start_time%month * 100 &
              + forcing_start_time%mday
        rtime = 0
        CALL taxisDefRdate(fdesc(ii)%out%taxisID,rdate)
        CALL taxisDefRtime(fdesc(ii)%out%taxisID,rtime)

        ! Assign the time axis to the variable list
        CALL vlistDefTaxis(fdesc(ii)%out%vlistID, fdesc(ii)%out%taxisID)

      ENDDO

    ENDIF
#endif
  END SUBROUTINE write_forcing_ini

!------------------------------------------------------------------------------

  SUBROUTINE write_forcing(fdesc, vdate, vtime)

    TYPE(field_desc), INTENT(INOUT) :: fdesc
    INTEGER,          INTENT(IN)    :: vdate, vtime
#ifndef NOCDI
    REAL(dp)                        :: forc_g(ie_g, je_g)

    INTEGER(i4)                     :: status, nmiss=0, year


    CHARACTER(LEN=maxnamelen) :: new_fname

    ! open file if necessary
    IF (p_pe == p_io) THEN
      year = vdate / 10000
      new_fname = spec_filename(fdesc%varname, year)  !  file name
      new_fname = TRIM(new_fname)//'_OUT'
      IF (fdesc%out%filename .ne. new_fname) THEN
        IF (TRIM(fdesc%out%filename) .ne. '') THEN
          CALL close_forcing_file(fdesc%out%filename,fdesc%out%streamID)
        END IF
        fdesc%out%filename = new_fname ! turn of the year is reached
        WRITE(0,*) 'Open output file ' // TRIM(new_fname)
        ! Create a dataset in netCDF format
        fdesc%out%streamID = streamOpenWrite(TRIM(new_fname), filetype_nc)
        IF ( fdesc%out%streamID < 0 )  &
          CALL stop_all ('Problem opening file '//TRIM(new_fname)//': '  &
               // cdiStringError(fdesc%out%streamID))
        ! Assign variable list to dataset
        CALL streamDefVlist(fdesc%out%streamID, fdesc%out%vlistID)
        fdesc%out%current_record = 0
      END IF
    ENDIF

    CALL gather(fdesc%dest,forc_g,p_io)

    IF (p_pe == p_io) THEN

      ! Set verification date
      CALL taxisDefVdate(fdesc%out%taxisID, vdate)

      ! Set verification time
      CALL taxisDefVtime(fdesc%out%taxisID, vtime)

      status = streamDefTimestep(fdesc%out%streamID, fdesc%out%current_record)

      CALL streamWriteVar(fdesc%out%streamID, fdesc%out%varID, &
           forc_g, nmiss)

      fdesc%out%current_record = fdesc%out%current_record + 1

    ENDIF

#endif
  END SUBROUTINE write_forcing

!------------------------------------------------------------------------------

  ELEMENTAL INTEGER FUNCTION record_idx(req_time)

    TYPE(time_desc), INTENT(IN) :: req_time

    record_idx = CEILING(REAL(model_seconds(req_time),dp) / forcing_frequency, i4)

  END FUNCTION record_idx

!------------------------------------------------------------------------------

  CHARACTER(LEN=maxnamelen) FUNCTION spec_filename(varname, year)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER,          INTENT(IN) :: year
    CHARACTER(LEN=4)             :: cfyear

    WRITE(cfyear,'(i4.4)') year

    spec_filename=TRIM(cforcdata)//'_'//TRIM(varname)//'_'//cfyear
    ! OMIP forcing data do not depend on year
    IF (TRIM(cforcdata) .EQ. 'OMIP' ) spec_filename=TRIM(cforcdata)//'_'//TRIM(varname)

  END FUNCTION spec_filename

!------------------------------------------------------------------------------

  INTEGER(i8) FUNCTION next_forcing_time_code(req_time, ref_year)

    TYPE(time_desc), INTENT(IN) :: req_time
    INTEGER, INTENT(IN) :: ref_year

    next_forcing_time_code = CEILING(REAL(model_seconds(req_time,ref_year),dp) &
                                          / forcing_frequency, i8)

  END FUNCTION next_forcing_time_code

!------------------------------------------------------------------------------

  SUBROUTINE time_interp_forcing(wgt1, wgt2, fdesc)
    REAL(dp), INTENT(IN)            :: wgt1, wgt2
    TYPE(field_desc), INTENT(INOUT) :: fdesc

    IF (wgt1 == 1._dp) THEN
      fdesc%dest = fdesc%ft1
    ELSEIF (wgt2 == 1._dp) THEN
      fdesc%dest = fdesc%ft2
    ELSE
      fdesc%dest = wgt1*fdesc%ft1 + wgt2*fdesc%ft2
    ENDIF

  END SUBROUTINE time_interp_forcing

!------------------------------------------------------------------------------

  SUBROUTINE runoff_ini(nlon, nlat, alat, alon, land)

    INTEGER(i4), INTENT(IN)    :: nlon, nlat
    REAL(dp),    INTENT(IN)    :: alat(nlon+2, nlat+2), alon(nlon+2,nlat+2)
    REAL(dp),    INTENT(INOUT) :: land(nlon+2, nlat+2)

    REAL(dp)    :: dist
    INTEGER(i4) :: i, j, ipos, jpos

    ALLOCATE(area(nlon,nlat), runoffij_cache(nlon,nlat))

    SELECT CASE (TRIM(cforcdata))
      CASE('OMIP') ! remove lakes
        land(41:58,111:124)=1._dp
    END SELECT

    WHERE (land .GE. 0.5_dp)
      land = 1._dp
    ELSEWHERE
      land = 0._dp
    ENDWHERE

    CALL surface_areas(nlon, nlat, alat(1,2:nlat+1), alon(2:nlon+1,1), area)

    DO i = 1, nlon
      DO j = 1, nlat
        CALL p_suchij(alat(i+1,j+1)*aradtogra, alon(i+1,j+1) * aradtogra, &
             1, ipos, jpos, dist, 1._wp, llocal_idx=.true.)

        runoffij_cache(i, j)%a = ipos
        runoffij_cache(i, j)%b = jpos
      END DO
    END DO

  END SUBROUTINE runoff_ini

!------------------------------------------------------------------------------

  SUBROUTINE runoff(nlon, nlat, ie, je, ri, riv, land)

    INTEGER(i4), INTENT(IN)    :: nlon, nlat, ie, je
    REAL(dp),    INTENT(INOUT) :: ri(nlon,nlat)
    REAL(dp),    INTENT(OUT)   :: riv(ie,je)
    REAL(dp),    INTENT(IN)    :: land(nlon,nlat)

    REAL(dp) :: rtest, rtest1

    INTEGER(i4) :: i, j, icount, ipos, jpos, jb

    rtest = 0._wp

    riv = 0._wp

    icount=0

    DO i=1, nlon
      DO j=1, nlat
        IF (land(i, j) .LT. 0.5_dp .AND. ri(i, j) .NE. 0._dp)THEN

          ri(i, j) = ri(i, j) * area(i, j)

          rtest = rtest + ri(i, j)

          icount = icount + 1

          ipos = runoffij_cache(i, j)%a
          jpos = runoffij_cache(i, j)%b
          IF (ipos /= -1) THEN
             jpos = runoffij_cache(i, j)%b
             riv(ipos, jpos) = riv(ipos, jpos) + ri(i, j)
          END IF
        END IF
      END DO
    END DO

    jb=MERGE(3, 2, lbounds_exch_tp .AND. p_joff == 0)
    rtest1 = SUM(riv(2:ie - 1, jb:je - 1))
    CALL global_sum(rtest1)

    IF (ldebug_forcing .AND. p_pe==p_io) &
         WRITE(0,'(a,2f20.10,i8)') 'runoff test: ', rtest, rtest1, icount

  END SUBROUTINE runoff

!------------------------------------------------------------------------------

  SUBROUTINE surface_areas(nlon, nlat, alat, alon, surfarea)

    ! Subroutine to calculate surface areas

    INTEGER(i4), INTENT(IN)  :: nlon, nlat
    REAL(dp),    INTENT(IN)  :: alon(nlon), alat(nlat)
    REAL(dp),    INTENT(OUT) :: surfarea(nlon, nlat)
    REAL(dp)                 :: alone, alonw, be, bn, bs, bw
    INTEGER(i4)              :: jx, jy

    DO jx = 1, nlon
      DO jy = 1, nlat

        bn = MERGE(api/2._dp, (alat(jy - 1) + alat(jy)) * 0.5_dp, jy .EQ. 1)
        bs = MERGE(-api/2._dp, (alat(jy + 1) + alat(jy)) * 0.5_dp, jy .EQ. nlat)
        alone = MERGE(alon(1), alon(jx + 1), jx .EQ. nlon)
        alonw = MERGE(alon(nlon), alon(jx - 1), jx .EQ. 1)

        IF (alone .LE. alon(jx)) alone = alone + 2._dp*api
        IF (alonw .GE. alon(jx)) alonw = alonw - 2._dp*api

        bw=MOD((alon(jx)+alonw)*0.5_dp, 2._dp*api)
        be=MOD((alon(jx)+alone)*0.5_dp, 2._dp*api)

        surfarea(jx, jy) = (be-bw) * radius  &
                         * (bn-bs) * radius  &
                         * COS(alat(jy))
      ENDDO
    ENDDO

  END SUBROUTINE surface_areas

!------------------------------------------------------------------------------

  SUBROUTINE destroy_objects(fdesc)
    TYPE(field_desc), INTENT(INOUT) :: fdesc
#ifndef NOCDI
    ! Destroy created objects
    CALL vlistDestroy(fdesc%out%vlistID)
    CALL taxisDestroy(fdesc%out%taxisID)
    !CALL gridDestroy (fdesc%out%gridID)
    !CALL zaxisDestroy(fdesc%out%zaxisID)
#endif
  END SUBROUTINE destroy_objects

!------------------------------------------------------------------------------

END MODULE mo_forcing
