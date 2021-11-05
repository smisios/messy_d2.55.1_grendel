! -*- f90 -*-
  ! -------------------------------------------------------------------
  !  CHANNEL INTERFACE FOR ECHAM5 / MPIOM OUTPUT
  ! -------------------------------------------------------------------

  !PRIVATE SUBROUTINE main_channel_initialize_gatts
  !PRIVATE SUBROUTINE main_channel_initialize_dims
  !PRIVATE SUBROUTINE main_channel_initialize_reprs
  !PRIVATE SUBROUTINE associate_streams_to_channels
  !PRIVATE SUBROUTINE reset_accu_stream_elements
  !PRIVATE SUBROUTINE timer_sync

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_initialize_gatts

    ! ECHAM5
    USE mo_advection,             ONLY: iadvec, no_advection, semi_lagrangian &
                                      , spitfire, tpcore

    ! MESSy
    USE messy_main_mpi_bi,          ONLY: p_bcast, p_io, p_parallel_io
    USE messy_main_data_bi,         ONLY: modstr_ECHAM5 => modstr &
                                        , modver_ECHAM5 => modver
    USE messy_main_grid_def_mem_bi, ONLY:  nn, nm, nk, lmidatm
    USE messy_main_timer,         ONLY: delta_time                     &
                                      , YEAR_START, MONTH_START, DAY_START &
                                      , HOUR_START, MINUTE_START, SECOND_START
    USE messy_main_channel,       ONLY: new_attribute, write_attribute &
                                      , AF_RST_CMP, EXP_NAME &
                                      , EXEC_CHECKSUM
    USE messy_main_compilerinfo_mem, ONLY: compiler_version, compiler_call &
                                         , compiler_flags &
                                         , compiler_cppdefs, compiler_includes
    USE messy_main_constants_mem,    ONLY: STRLEN_ULONG

    IMPLICIT NONE

    !  EXTERNAL SUBROUTIONES
    EXTERNAL util_os_system, util_user_name, util_node_name

    INTRINSIC :: DATE_AND_TIME, TRIM, LEN_TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_initialize_gatts'
    INTEGER :: status
    INTEGER :: nlena, nlenb, nlenc
    CHARACTER (256) :: os_name, user_name, host_name
    CHARACTER (8)   :: ydate
    CHARACTER (10)  :: ytime
    CHARACTER (DATE_TIME_STR_LEN) :: str_date_time
    CHARACTER (LEN=2)            :: istr
    INTEGER :: i1, i2, il, ic

    ! INIT
    IF (p_parallel_io) THEN
       CALL util_os_system (os_name, nlena)
       CALL util_user_name (user_name, nlenb)
       CALL util_node_name (host_name, nlenc)
    END IF
    CALL p_bcast (os_name, p_io)
    CALL p_bcast (nlena, p_io)
    CALL p_bcast (user_name, p_io)
    CALL p_bcast (nlenb, p_io)
    CALL p_bcast (host_name, p_io)
    CALL p_bcast (nlenc, p_io)

    ! avoid deviating entries between different PEs;
    ! pnetcdf gives an error message about deviating
    ! attributes
    IF (p_parallel_io) THEN
     CALL DATE_AND_TIME(ydate, ytime)
    END IF
    CALL p_bcast (ydate, p_io)
    CALL p_bcast (ytime, p_io)

    ! EXPERIMENT NAME
    CALL new_attribute(status, 'MESSy_experiment', c=TRIM(EXP_NAME))
    CALL channel_halt(substr, status)

    ! CHECKSUM OF EXECUTABLE
    CALL new_attribute(status, 'EXEC_CHECKSUM', c=TRIM(EXEC_CHECKSUM))
    CALL channel_halt(substr, status)

    ! ECHAM5 SETUP
    ! - VERSION
    CALL new_attribute(status, 'GCM' &
         , c = modstr_ECHAM5//' version '//modver_ECHAM5//&
         &', Max-Planck Institute for Meteorology, Hamburg' )
    CALL channel_halt(substr, status)
    ! - TRUNCATION
    CALL new_attribute(status, 'GCM_spherical_trunc_n', i=nn, iflag=AF_RST_CMP)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'GCM_spherical_trunc_m', i=nm, iflag=AF_RST_CMP)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'GCM_spherical_trunc_k', i=nk, iflag=AF_RST_CMP)
    CALL channel_halt(substr, status)
    ! - MODI
    IF (lmidatm) THEN
       CALL new_attribute(status, 'GCM_vertical_mode' &
            , c='middle atmosphere (MA)', iflag=AF_RST_CMP)
    ELSE
       CALL new_attribute(status, 'GCM_vertical_mode' &
            , c='normal', iflag=AF_RST_CMP)
    END IF
    CALL channel_halt(substr, status)
    !
    CALL new_attribute(status, 'GCM_horizontal_mode', c='global' &
         , iflag=AF_RST_CMP)
    CALL channel_halt(substr, status)
    ! - ADVECTION
    SELECT CASE(iadvec)
    CASE(no_advection)
       CALL new_attribute(status, 'GCM_advection', c='no advection')
    CASE(semi_lagrangian)
       CALL new_attribute(status, 'GCM_advection', c='Semi-Lagrangian')
    CASE(spitfire)
       CALL new_attribute(status, 'GCM_advection', c='SPITFIRE')
    CASE(tpcore)
       CALL new_attribute(status, 'GCM_advection', c='Lin&Rood')
    CASE DEFAULT
       CALL new_attribute(status, 'GCM_advection', c='unknown')
    END SELECT
    CALL channel_halt(substr, status)
    ! - SIMULATION TIME
    WRITE(str_date_time,'(i4,i2.2,i2.2,a1,3(i2.2))') &
         YEAR_START, MONTH_START, DAY_START,         &
         ' ',HOUR_START, MINUTE_START, SECOND_START
    CALL new_attribute(status, 'GCM_start_date_time' &
         , c=str_date_time)
    CALL channel_halt(substr, status)
    ! - TIMESTEP LENGTH
    CALL new_attribute(status, 'GCM_timestep' &
         , r=delta_time)
    CALL channel_halt(substr, status)

    ! COMPILER INFORMATION
    il = LEN_TRIM(compiler_version)
    IF (il > STRLEN_ULONG) THEN
       i1 = 1
       i2 = STRLEN_ULONG
       ic = 0
       DO
          ic = ic + 1
          WRITE(istr,'(I2.2)') ic
          CALL new_attribute(status, 'F95_COMPILER_VERSION_'//istr &
               , c=compiler_version(i1:i2))
          CALL channel_halt(substr, status)
          IF (i2 == il) EXIT
          i1 = i2 + 1
          i2 = MIN(i2+STRLEN_ULONG, il)
       END DO
    ELSE
       CALL new_attribute(status, 'F95_COMPILER_VERSION' &
            , c=TRIM(compiler_version))
       CALL channel_halt(substr, status)
    ENDIF

    il = LEN_TRIM(compiler_call)
    IF (il > STRLEN_ULONG) THEN
       i1 = 1
       i2 = STRLEN_ULONG
       ic = 0
       DO
          ic = ic + 1
          WRITE(istr,'(I2.2)') ic
          CALL new_attribute(status, 'F95_COMPILER_CALL_'//istr &
               , c=compiler_call(i1:i2))
          CALL channel_halt(substr, status)
          IF (i2 == il) EXIT
          i1 = i2 + 1
          i2 = MIN(i2+STRLEN_ULONG, il)
       END DO
    ELSE
       CALL new_attribute(status, 'F95_COMPILER_CALL', c=TRIM(compiler_call))
       CALL channel_halt(substr, status)
    END IF

    il = LEN_TRIM(compiler_flags)
    IF (il > STRLEN_ULONG) THEN
       i1 = 1
       i2 = STRLEN_ULONG
       ic = 0
       DO
          ic = ic + 1
          WRITE(istr,'(I2.2)') ic
          CALL new_attribute(status, 'F95_COMPILER_FLAGS_'//istr &
               , c=compiler_flags(i1:i2))
          CALL channel_halt(substr, status)
          IF (i2 == il) EXIT
          i1 = i2 + 1
          i2 = MIN(i2+STRLEN_ULONG, il)
       END DO
    ELSE
       CALL new_attribute(status, 'F95_COMPILER_FLAGS', c=TRIM(compiler_flags))
       CALL channel_halt(substr, status)
    END IF

    il = LEN_TRIM(compiler_cppdefs)
    IF (il > STRLEN_ULONG) THEN
       i1 = 1
       i2 = STRLEN_ULONG
       ic = 0
       DO
          ic = ic + 1
          WRITE(istr,'(I2.2)') ic
          CALL new_attribute(status, 'F95_PREPROC_DEFINITIONS'//istr &
               , c=compiler_cppdefs(i1:i2))
          CALL channel_halt(substr, status)
          IF (i2 == il) EXIT
          i1 = i2 + 1
          i2 = MIN(i2+STRLEN_ULONG, il)
       END DO
    ELSE
       CALL new_attribute(status, 'F95_PREPROC_DEFINITIONS' &
            , c=TRIM(compiler_cppdefs))
       CALL channel_halt(substr, status)
    END IF

    il = LEN_TRIM(compiler_includes)
    IF (il > STRLEN_ULONG) THEN
       i1 = 1
       i2 = STRLEN_ULONG
       ic = 0
       DO
          ic = ic + 1
          WRITE(istr,'(I2.2)') ic
          CALL new_attribute(status, 'F95_COMPILER_INCLUDES_'//istr &
               , c=compiler_includes(i1:i2))
          CALL channel_halt(substr, status)
          IF (i2 == il) EXIT
          i1 = i2 + 1
          i2 = MIN(i2+STRLEN_ULONG, il)
       END DO
    ELSE
       CALL new_attribute(status, 'F95_COMPILER_INCLUDES' &
            , c=TRIM(compiler_includes))
       CALL channel_halt(substr, status)
    END IF

    ! OPERATION
    CALL new_attribute(status, 'operating_date_time' &
         , c = ydate(1:8)//' '//ytime(1:6) )
    CALL channel_halt(substr, status)

    CALL new_attribute(status, 'operating_system', c=os_name(1:nlena) )
    CALL channel_halt(substr, status)

    CALL new_attribute(status, 'operating_host', c=host_name(1:nlenc) )
    CALL channel_halt(substr, status)

    CALL new_attribute(status, 'operating_user', c=user_name(1:nlenb) )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) THEN
       CALL write_attribute(status)
       CALL channel_halt(substr, status)
    END IF

  END SUBROUTINE main_channel_initialize_gatts
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_initialize_dims

    ! MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_grid_def_mem_bi,  ONLY: nvclev, vct    &  ! hyai, hybi
                                         , nlev           &  ! number of levels
                                         , nlon, ngl      &  !
                                         , nsp            &
                                         , nhgl, nmp1
    USE messy_main_grid_def_bi,      ONLY: philon, philat
    ! MESSy/SMCL
    USE messy_main_timer,     ONLY: &
                                    delta_time     &
                                  , current_time_step &
                                  , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND &
                                  , YEAR_START, MONTH_START, DAY_START     &
                                  , HOUR_START, MINUTE_START, SECOND_START
    USE messy_main_channel_dimensions,   ONLY: new_dimension            &
                                             , write_dimension          &
                                             , add_dimension_variable   &
                                             , add_dimension_variable_att

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_initialize_dims'
    INTEGER                             :: status
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: array
    INTEGER                             :: i, n
    CHARACTER(LEN=80)                   :: start_date_string
    REAL(dp)                            :: yyyymmdd

    ! #############
    ! ### TIME ####
    ! #############
    CALL new_dimension(status, DIMID_TIME, 'time', 1, .TRUE.)
    CALL channel_halt(substr, status)
    ! - TIME
    CALL add_dimension_variable(status, 'time', 'time', (/ 0.0_DP /))
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'time', 'time', &
         'long_name', c='time')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'time', 'time', &
         'bounds', c='time_bnds')
    CALL channel_halt(substr, status)
    !
    WRITE(start_date_string, &
         '("day since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
         YEAR_START,MONTH_START,DAY_START,HOUR_START,MINUTE_START,SECOND_START
    CALL add_dimension_variable_att(status, 'time', 'time', &
         'units', c=start_date_string)
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'time', 'time', &
         'calendar', c='gregorian')
    CALL channel_halt(substr, status)
    ! - YYYYMMDD
    yyyymmdd = ABS(YEAR)*10000 + MONTH*100 + DAY  &
         + (HOUR*3600 + MINUTE*60 + SECOND)/86400.
    IF (YEAR<0) yyyymmdd = -yyyymmdd
    CALL add_dimension_variable(status, 'time', 'YYYYMMDD', (/ yyyymmdd /))
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'time', 'YYYYMMDD', &
         'long_name', c='time')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'time', 'YYYYMMDD', &
         'units', c='days as %Y%m%d.%f')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'time', 'YYYYMMDD', &
         'calendar', c='gregorian')
    CALL channel_halt(substr, status)
    ! - DT
    CALL add_dimension_variable(status, 'time', 'dt', (/ delta_time /))
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'time', 'dt', &
         'long_name', c='delta_time')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'time', 'dt', &
         'units', c='s')
    CALL channel_halt(substr, status)
    ! - CURRENT TIME STEP
    CALL add_dimension_variable(status, 'time', 'nstep' &
         , (/ REAL(current_time_step, DP) /))
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'time', 'nstep', &
         'long_name', c='current time step')
    CALL channel_halt(substr, status)

    ! ##################
    ! ### GRIDPOINT ####
    ! ##################

    ! MID-LAYERS (FULL LEVELS)
    ALLOCATE(array(nlev))
    !
    CALL new_dimension(status, DIMID_LEV, 'lev', nlev)
    CALL channel_halt(substr, status)
    !
    DO i=1, nlev
       array(i) = REAL(i, DP)
    END DO
    CALL add_dimension_variable(status, 'lev', 'lev', array)
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'lev', 'lev', &
         'long_name', c='hybrid level at layer midpoints')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'lev', 'lev', &
         'standard_name', c='hybrid_sigma_pressure')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'lev', 'lev', &
         'units', c='level')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'lev', 'lev', &
         'positive', c='down')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'lev', 'lev', &
         'formula', c='hyam hybm (press=hyam+hybm*aps)')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'lev', 'lev', &
         'borders', c='ilev')
    CALL channel_halt(substr, status)

    array(:) = (vct(1:nvclev-1) + vct(2:nvclev))/2
    CALL add_dimension_variable(status, 'lev', 'hyam', array)
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'lev', 'hyam', &
         'long_name', c='hybrid A coefficient at layer midpoints')
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'lev', 'hyam', &
         'units', c='Pa')
    CALL channel_halt(substr, status)

    array(:) = (vct(nvclev+1:2*nvclev-1) + vct(nvclev+2:2*nvclev))/2
    CALL add_dimension_variable(status, 'lev', 'hybm', array)
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'lev', 'hybm', &
         'long_name', c='hybrid B coefficient at layer midpoints')
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'lev', 'hybm', &
         'units', c='1')
    CALL channel_halt(substr, status)

    DEALLOCATE(array)

    ! INTERFCAE-LAYERS (HALF LEVELS)
    ALLOCATE(array(nlev+1))
    !
    CALL new_dimension(status, DIMID_ILEV, 'ilev', nlev+1)
    CALL channel_halt(substr, status)
    !
    DO i=1, nlev+1
       array(i) = REAL(i, DP)
    END DO
    CALL add_dimension_variable(status, 'ilev', 'ilev', array)
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'ilev', 'ilev', &
         'long_name', c='hybrid level at layer interfaces')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'ilev', 'ilev', &
         'standard_name', c='hybrid_sigma_pressure')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'ilev', 'ilev', &
         'units', c='level')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'ilev', 'ilev', &
         'positive', c='down')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'ilev', 'ilev', &
         'formula', c='hyai hybi (ilev=hyai+hybi*aps)')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'ilev', 'ilev', &
         'formula_terms', c='ap: hyai b: hybi ps: aps')
    CALL channel_halt(substr, status)

    array(:) = vct(1:nvclev)
    CALL add_dimension_variable(status, 'ilev', 'hyai', array)
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'ilev', 'hyai', &
         'long_name', c='hybrid A coefficient at layer interfaces')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'ilev', 'hyai', &
         'units', c='Pa')
    CALL channel_halt(substr, status)

    array(:) = vct(nvclev+1:2*nvclev)
    CALL add_dimension_variable(status, 'ilev', 'hybi', array)
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'ilev', 'hybi', &
         'long_name', c='hybrid B coefficient at layer interfaces')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'ilev', 'hybi', &
         'units', c='1')
    CALL channel_halt(substr, status)

    DEALLOCATE(array)

    ! LONGITUDES
    n = nlon
    ALLOCATE(array(n))
    array(:) = philon(1:nlon)

    CALL new_dimension(status, DIMID_LON, 'lon', n)
    CALL channel_halt(substr, status)
    CALL add_dimension_variable(status, 'lon', 'lon', array)
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'lon', 'lon', &
         'long_name', c='longitude')
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'lon', 'lon', &
         'units', c='degrees_east')
    CALL channel_halt(substr, status)

    DEALLOCATE(array)

    ! LATITUDES
    n = ngl
    ALLOCATE(array(n))
    array(:) = philat(1:ngl)

    CALL new_dimension(status, DIMID_LAT, 'lat', n)
    CALL channel_halt(substr, status)
    CALL add_dimension_variable(status, 'lat', 'lat', array)
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'lat', 'lat', &
         'long_name', c='latitude')
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'lat', 'lat', &
         'units', c='degrees_north')
    CALL channel_halt(substr, status)

    DEALLOCATE(array)

    ! #################
    ! ### SPECTRAL ####
    ! #################

    ! COMPLEX
    CALL new_dimension(status, DIMID_COMPLEX, 'complex', 2)
    CALL channel_halt(substr, status)
    !
    ! SPECTRAL COEFFICIENTS
    CALL new_dimension(status, DIMID_NSP, 'nsp', nsp)
    CALL channel_halt(substr, status)

    ! ############################################
    ! ### FOURIER (ANTI-SYMMETRIC, SYMMETRIC) ####
    ! ############################################

    ! HALF NUMBER OF LATITUDES
    CALL new_dimension(status, DIMID_NHGL, 'nhgl', nhgl)
    CALL channel_halt(substr, status)

    ! ZONAL WAVE NUMBER + 1
    CALL new_dimension(status, DIMID_NMP1, 'nmp1', nmp1)
    CALL channel_halt(substr, status)

    ! ################
    ! ### SPECIAL ####
    ! ################

    ! BELOWSURFACE
    CALL new_dimension(status, DIMID_BELOWSF, 'belowsf', 5)
    CALL channel_halt(substr, status)

    CALL add_dimension_variable(status, 'belowsf', 'belowsf' &
         , (/0.03_dp, 0.19_dp, 0.78_dp, 2.68_dp, 6.98_dp/) )
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'belowsf', 'belowsf', &
         'long_name', c='levels below surface')
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'belowsf', 'belowsf', &
         'units', c='m')
    CALL channel_halt(substr, status)

    ! 2 LEVELS
    CALL new_dimension(status, DIMID_2LEV, 'l2', 2)
    CALL channel_halt(substr, status)

    ! 2 boundaries
     CALL new_dimension(status, DIMID_TBNDS, 'tbnds', tbnds)
     CALL channel_halt(substr, status)

    ! #################
    ! DIAGNOSTIC OUTPUT
    ! #################
    IF (p_parallel_io) THEN
       CALL write_dimension(status)
       CALL channel_halt(substr, status)
    END IF

  END SUBROUTINE main_channel_initialize_dims
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_initialize_reprs

    ! MESSy
    USE messy_main_mpi_bi,               ONLY: p_parallel_io, p_pe, p_nprocs &
                                             , dcl, dcg
    USE messy_main_grid_def_mem_bi,      ONLY: nproma, ngpblks         &
                                             , snsp, nllev, nllevp1    &
                                             , nlm
    USE messy_main_channel_repr,         ONLY: new_representation        &
                                             , write_representation_dc   &
                                             , set_representation_decomp &
                                             , write_representation      &
                                             , AUTO, IRANK               &
                                             , PIOTYPE_COL

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_initialize_reprs'
    INTEGER                     :: status
    INTEGER, DIMENSION(:), POINTER :: scdim => NULL()
    !
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg    = 0  ! max no. of segments
    INTEGER                          :: nseg_pe = 0  ! segments on this PE
    INTEGER                          :: n            ! counter
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()
    INTEGER                          :: mp, mp1, im, np, mpgl ! for SPECTRAL
    INTEGER                          :: ke, nk ! for FAS

    ! ##################
    ! ### GRIDPOINT ####
    ! ##################

    nseg = gp_nseg
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    start(:,:) = 1
    cnt(:,:) = 1
    meml(:,:) = 1
    memu(:,:) = 1

    ! ---------------------------------------------

    CALL new_representation(status, GP_3D_MID, 'GP_3D_MID' &
         , rank = 3, link = 'xxx-', dctype = DC_GP               &
         , dimension_ids = (/ DIMID_LON, DIMID_LEV, DIMID_LAT /) &
         , ldimlen       = (/ nproma   , AUTO     , ngpblks   /) &
         , output_order  = (/ 1,3,2 /)                           &
         , axis = 'XZY-'                                         &
         )
    CALL channel_halt(substr, status)

    start(:,1) = dcl%glons(:)
    start(:,2) = 1
    start(:,3) = dcl%glats(:)

    cnt(:,1) = dcl%glone(:)-dcl%glons(:)+1
    cnt(:,2) = dcl%nlev
    cnt(:,3) = dcl%glate(:)-dcl%glats(:)+1

    meml(:,1) = 1
    meml(:,2) = 1
    meml(:,3) = (/ 1, dcl%nglh(1)+1 /)

    memu(:,1) = dcl%nglon
    memu(:,2) = dcl%nlev
    memu(:,3) = (/ dcl%nglh(1), dcl%nglh(1)+dcl%nglh(2)/)

    CALL set_representation_decomp(status, GP_3D_MID &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! SAVE GP DECOMPOSITION INFORAMTION FOR USAGE IN SUBMODELS
    nseg = gp_nseg
    ALLOCATE(gp_start(nseg,IRANK))
    ALLOCATE(gp_cnt(nseg,IRANK))
    ALLOCATE(gp_meml(nseg,IRANK))
    ALLOCATE(gp_memu(nseg,IRANK))

    gp_start(:,1:2) = start(:,1:2)
    gp_cnt(:,1:2) = cnt(:,1:2)
    gp_meml(:,1:2) = meml(:,1:2)
    gp_memu(:,1:2) = memu(:,1:2)

    gp_start(:,4) = start(:,3)
    gp_cnt(:,4) = cnt(:,3)
    gp_meml(:,4) = meml(:,3)
    gp_memu(:,4) = memu(:,3)

    gp_start(:,3) = 1
    gp_cnt(:,3) = 1
    gp_meml(:,3) = 1
    gp_memu(:,3) = 1

    ! -----------------------------------------------------

    CALL new_representation(status, GP_3D_INT, 'GP_3D_INT' &
         , rank = 3, link = 'xxx-', dctype = DC_GP                &
         , dimension_ids = (/ DIMID_LON, DIMID_ILEV, DIMID_LAT /) &
         , ldimlen       = (/ nproma   , AUTO      , ngpblks   /) &
         , output_order  = (/ 1,3,2 /)                            &
         , axis = 'XZY-'                                          &
         )
    CALL channel_halt(substr, status)

    cnt(:,2)  = dcl%nlev+1
    memu(:,2) = dcl%nlev+1

    CALL set_representation_decomp(status, GP_3D_INT &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -----------------------------------------------------

    CALL new_representation(status, GP_2D_HORIZONTAL, 'GP_2D_HORIZONTAL' &
         , rank = 2, link = 'xx--', dctype = DC_GP       &
         , dimension_ids = (/ DIMID_LON, DIMID_LAT /)    &
         , ldimlen       = (/ nproma   , ngpblks   /)    &
         , axis = 'XY--'                                 &
         )
    CALL channel_halt(substr, status)

    start(:,2) = start(:,3)
    cnt(:,2)   = cnt(:,3)
    meml(:,2)  = meml(:,3)
    memu(:,2)  = memu(:,3)

    start(:,3) = 1
    cnt(:,3)   = 1
    meml(:,3)  = 1
    memu(:,3)  = 1

    CALL set_representation_decomp(status, GP_2D_HORIZONTAL &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -----------------------------------------------------

    CALL new_representation(status, GP_3D_1LEV, 'GP_3D_1LEV' &
         , rank = 2, link = 'x-x-', dctype = DC_GP       &
         , dimension_ids = (/ DIMID_LON, DIMID_LAT /)    &
         , ldimlen       = (/ nproma   , ngpblks   /)    &
         , axis = 'X-Y-'                                 &
         )
    CALL channel_halt(substr, status)

    ! NOTE: THE RANK IS W.R.T. MEMORY, I.E., ACCORDING TO 'link' !!!

    start(:,3) = start(:,2)
    cnt(:,3)   = cnt(:,2)
    meml(:,3)  = meml(:,2)
    memu(:,3)  = memu(:,2)

    start(:,2) = 1
    cnt(:,2)   = 1
    meml(:,2)  = 1
    memu(:,2)  = 1

    CALL set_representation_decomp(status, GP_3D_1LEV &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -----------------------------------------------------

    CALL new_representation(status, GP_1D_LEV, 'GP_1D_LEV' &
         , rank = 1, link = '-x--', dctype = DC_BC         &
         , dimension_ids = (/ DIMID_LEV /)                 &
         , ldimlen       = (/ AUTO   /)                    &
         , axis = '-Z--'                                   &
         )
    CALL channel_halt(substr, status)

    ! NOTE: THE RANK IS W.R.T. MEMORY, I.E., ACCORDING TO 'link' !!!

    start(:,:) = 1
    cnt(:,:)   = 1
    meml(:,:)  = 1
    memu(:,:)  = 1

    start(:,2) = 1
    cnt(:,2)   = dcl%nlev
    meml(:,2)  = 1
    memu(:,2)  = dcl%nlev

    CALL set_representation_decomp(status, GP_1D_LEV &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -----------------------------------------------------

    CALL new_representation(status, GP_1D_ILEV, 'GP_1D_ILEV' &
         , rank = 1, link = '-x--', dctype = DC_BC           &
         , dimension_ids = (/ DIMID_ILEV /)                  &
         , ldimlen       = (/ AUTO   /)                      &
         , axis = '-Z--'                                     &
         )
    CALL channel_halt(substr, status)

    ! NOTE: THE RANK IS W.R.T. MEMORY, I.E., ACCORDING TO 'link' !!!

    start(:,:) = 1
    cnt(:,:)   = 1
    meml(:,:)  = 1
    memu(:,:)  = 1

    start(:,2) = 1
    cnt(:,2)   = dcl%nlev +1
    meml(:,2)  = 1
    memu(:,2)  = dcl%nlev +1

    CALL set_representation_decomp(status, GP_1D_ILEV &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -----------------------------------------------------

    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    ! -----------------------------------------------------

    ! #################
    ! ### SPECTRAL ####
    ! #################

    ! MAX NUMBER OF SEGMENTS
    nseg = 0
    DO n=1, p_nprocs
       IF (dcg(n)%nsm > nseg) nseg = dcg(n)%nsm
    END DO

    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    nseg_pe = dcl%nsm

    DO im=1,nseg
       start(im,1) = 1
       cnt(im,1)   = 0
       meml(im,1)  = 1
       memu(im,1)  = 1

       start(im,2) = 1
       cnt(im,2)   = 2
       meml(im,2)  = 1
       memu(im,2)  = 2

       start(im,3) = 1
       cnt(im,3)   = 0
       meml(im,3)  = 1
       memu(im,3)  = 0

       start(im,4) = 0
       cnt(im,4)   = 0
       meml(im,4)  = 0
       memu(im,4)  = 0
    END DO

    ! -----------------------------------------------------

    CALL new_representation(status, SP_3D_MID, 'SP_3D_MID' &
         , rank = 3, link = 'xxx-', dctype = DC_SP         &
         , dimension_ids = (/ DIMID_LEV, DIMID_COMPLEX, DIMID_NSP /) &
         , ldimlen       = (/ AUTO     , AUTO         , snsp      /) &
         , output_order  = (/ 2,3,1 /)                               &
         , axis = 'Z---'                                             &
         )
    CALL channel_halt(substr, status)

    mp=0
    DO im=1,nseg_pe
       mp1  = dcl%sm(im)+1
       np   = dcl%snnp(im)
       mpgl = dcl%nmp(mp1) + dcl%snn0(im)

       start(im,1) = 1
       cnt(im,1)   = dcl%nlev
       meml(im,1)  = 1
       memu(im,1)  = dcl%nlev

       start(im,2) = 1
       cnt(im,2)   = 2
       meml(im,2)  = 1
       memu(im,2)  = 2

       start(im,3) = mpgl+1
       cnt(im,3)   = np
       meml(im,3)  = mp+1
       memu(im,3)  = mp+np

       mp = mp + np
    END DO

    CALL set_representation_decomp(status, SP_3D_MID &
         , start, cnt, memu, meml, piotype=PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -----------------------------------------------------

    CALL new_representation(status, SP_3D_INT, 'SP_3D_INT' &
         , rank = 3, link = 'xxx-', dctype = DC_SP         &
         , dimension_ids = (/ DIMID_ILEV, DIMID_COMPLEX, DIMID_NSP /) &
         , ldimlen       = (/ AUTO      , AUTO         , snsp      /) &
         , output_order  = (/ 2,3,1 /)                                &
         , axis = 'Z---'                                              &
         )
    CALL channel_halt(substr, status)

    DO im=1,nseg_pe
       cnt(im,1)   = dcl%nlev+1
       memu(im,1)  = dcl%nlev+1
    END DO

    CALL set_representation_decomp(status, SP_3D_INT &
         , start, cnt, memu, meml, piotype=PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -----------------------------------------------------

    CALL new_representation(status, SP_2D_HORIZONTAL, 'SP_2D_HORIZONTAL' &
         , rank = 2, link = '-xx-', dctype = DC_SP        &
         , dimension_ids = (/ DIMID_COMPLEX, DIMID_NSP /) &
         , ldimlen       = (/ AUTO         , snsp      /) &
         !, output_order  = (/ 2,3 /)                      &
         , axis = '----'                                  &
         )
    CALL channel_halt(substr, status)

    ! NOTE: THE RANK IS W.R.T. MEMORY, I.E., ACCORDING TO 'link' !!!

    DO im=1,nseg_pe
       start(im,1) = 1
       cnt(im,1)   = 1
       meml(im,1)  = 1
       memu(im,1)  = 1
    END DO

    CALL set_representation_decomp(status, SP_2D_HORIZONTAL &
         , start, cnt, memu, meml, piotype=PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -----------------------------------------------------

    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    ! -----------------------------------------------------

    ! ############################################
    ! ### FOURIER (ANTI-SYMMETRIC, SYMMETRIC) ####
    ! ############################################

    ! -----------------------------------------------------

    ! MAX NUMBER OF SEGMENTS
    nseg = 0
    DO n=1, p_nprocs
       IF (dcg(n)%nlm > nseg) nseg = dcg(n)%nlm
    END DO

    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    nseg_pe = dcl%nlm
    start(:,:) = 1
    cnt(:,:) = 0
    meml(:,:) = 1
    memu(:,:) = 1

    ! -----------------------------------------------------

    CALL new_representation(status, FAS_MID, 'FAS_MID' &
         , rank = 4, link = 'xxxx', dctype = DC_SA        &
         , dimension_ids = (/ DIMID_LEV,   DIMID_COMPLEX  &
                            , DIMID_NMP1,  DIMID_NHGL /)  &
         , ldimlen       = (/ nllev, AUTO, nlm , AUTO /)  &
         , axis = 'Z---'                                  &
         )
    CALL channel_halt(substr, status)

    ke = MIN (dcl%lleve,dcl%nlev)
    nk = MAX(ke - dcl%llevs + 1, 0)

    DO im=1,nseg_pe
       mp1  = dcl%lm(im)+1

       start(im,1) = dcl%llevs
       cnt(im,1)   = nk
       meml(im,1)  = 1
       memu(im,1)  = nk

       start(im,2) = 1
       cnt(im,2)   = 2
       meml(im,2)  = 1
       memu(im,2)  = 2

       start(im,3) = mp1
       cnt(im,3)   = 1
       meml(im,3)  = im
       memu(im,3)  = im

       start(im,4) = 1
       cnt(im,4)   = dcl%nlat/2
       meml(im,4)  = 1
       memu(im,4)  = dcl%nlat/2
    END DO

    CALL set_representation_decomp(status, FAS_MID &
         , start, cnt, memu, meml, piotype=PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -----------------------------------------------------

    CALL new_representation(status, FAS_INT, 'FAS_INT' &
         , rank = 4, link = 'xxxx', dctype = DC_SA          &
         , dimension_ids = (/ DIMID_ILEV,  DIMID_COMPLEX    &
                            , DIMID_NMP1,  DIMID_NHGL   /)  &
         , ldimlen       = (/ nllevp1, AUTO, nlm , AUTO /)  &
         , axis = 'Z---'                                    &
         )
    CALL channel_halt(substr, status)

    ke = MIN (dcl%lleve,dcl%nlev+1)
    nk = MAX(ke - dcl%llevs + 1, 0)

    DO im=1,nseg_pe
       start(im,1) = dcl%llevs
       cnt(im,1)   = nk
       meml(im,1)  = 1
       memu(im,1)  = nk
    END DO

    CALL set_representation_decomp(status, FAS_INT &
         , start, cnt, memu, meml, piotype=PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -----------------------------------------------------

    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    ! -----------------------------------------------------

    nseg = 1
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    start(:,:) = 1
    cnt(:,:) = 1
    meml(:,:) = 1
    memu(:,:) = 1

    ! -----------------------------------------------------

    CALL new_representation(status, FAS_MID_ZM, 'FAS_MID_ZM' &
         , rank = 2, link = 'xx--', dctype = DC_SA        &
         , dimension_ids = (/ DIMID_LEV, DIMID_NHGL   /)  &
         , ldimlen       = (/ nllev, AUTO             /)  &
         , axis = 'Z---'                                  &
         )
    CALL channel_halt(substr, status)

    ke = MIN (dcl%lleve,dcl%nlev)
    nk = MAX(ke - dcl%llevs + 1, 0)

    start(:,1) = dcl%llevs
    cnt(:,1)   = nk
    meml(:,1)  = 1
    memu(:,1)  = nk

    start(:,2) = 1
    cnt(:,2)   = dcl%nlat/2
    meml(:,2)  = 1
    memu(:,2)  = dcl%nlat/2

    CALL set_representation_decomp(status, FAS_MID_ZM &
         , start, cnt, memu, meml, piotype=PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -----------------------------------------------------

    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    ! -----------------------------------------------------

    ! ###############
    ! ### SPECIAL ###
    ! ###############

    nseg = 2
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    start(:,:) = 1
    cnt(:,:) = 1
    meml(:,:) = 1
    memu(:,:) = 1

    ! -----------------------------------------------------

    CALL new_representation(status, GP_3D_BELOWSF, 'GP_3D_BELOWSF'   &
         , rank = 3, link = 'xxx-', dctype = DC_GP                   &
         , dimension_ids = (/ DIMID_LON, DIMID_BELOWSF, DIMID_LAT /) &
         , ldimlen       = (/ nproma   , AUTO     , ngpblks   /)     &
         , output_order  = (/ 1,3,2 /)                               &
         , axis = 'XZY-'                                             &
         )
    CALL channel_halt(substr, status)

    start(:,1) = dcl%glons(:)
    start(:,2) = 1
    start(:,3) = dcl%glats(:)

    cnt(:,1) = dcl%glone(:)-dcl%glons(:)+1
    cnt(:,2) = 5
    cnt(:,3) = dcl%glate(:)-dcl%glats(:)+1

    meml(:,1) = 1
    meml(:,2) = 1
    meml(:,3) = (/ 1, dcl%nglh(1)+1 /)

    memu(:,1) = dcl%nglon
    memu(:,2) = 5
    memu(:,3) = (/ dcl%nglh(1), dcl%nglh(1)+dcl%nglh(2)/)

    CALL set_representation_decomp(status, GP_3D_BELOWSF &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -----------------------------------------------------

    CALL new_representation(status, GP_3D_2LEV, 'GP_3D_2LEV'         &
         , rank = 3, link = 'xxx-', dctype = DC_GP                   &
         , dimension_ids = (/ DIMID_LON, DIMID_2LEV, DIMID_LAT /)    &
         , ldimlen       = (/ nproma   , AUTO     , ngpblks   /)     &
         , output_order  = (/ 1,3,2 /)                               &
         , axis = 'XZY-'                                             &
         )
    CALL channel_halt(substr, status)

    cnt(:,2)  = 2
    memu(:,2) = 2

    CALL set_representation_decomp(status, GP_3D_2LEV &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -----------------------------------------------------

    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    ! -----------------------------------------------------

    ALLOCATE(scdim(0))
    scdim = 0
    CALL new_representation(status, SCALAR, 'SCALAR'                 &
         , rank = 0, link = '----', dctype = DC_BC                   &
         , dimension_ids = scdim &
         , ldimlen       = scdim &
         , axis = '----'                                             &
         )
    CALL channel_halt(substr, status)
    DEALLOCATE(scdim)
    NULLIFY(scdim)

    nseg = 1
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    start(:,:) = 1
    cnt(:,:) = 1
    meml(:,:) = 1
    memu(:,:) = 1

    CALL set_representation_decomp(status, SCALAR &
         , start, cnt, memu, meml, piotype=PIOTYPE_COL)
    CALL channel_halt(substr, status)

    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

   ! -----------------------------------------------------

    CALL new_representation(status, REPR_TIMEBNDS, 'TIME_BNDS' &
         , rank = 1, link = 'x---', dctype = DC_BC             &
         , dimension_ids = (/DIMID_TBNDS/)                     &
         , ldimlen       = (/ tbnds     /)                     &
         , axis = 'N---'                                       &
         )
    CALL channel_halt(substr, status)

    nseg = 1
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    start(:,:) = 1
    cnt(:,:)   = 1
    meml(:,:)  = 1
    memu(:,:)  = 1

    cnt(:,1)  = tbnds
    memu(:,1) = tbnds

    CALL set_representation_decomp(status, REPR_TIMEBNDS &
         , start, cnt, memu, meml, piotype=PIOTYPE_COL)
    CALL channel_halt(substr, status)

    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

   ! -----------------------------------------------------

    ! #################
    ! DIAGNOSTIC OUTPUT
    ! #################
    IF (p_parallel_io) THEN
       CALL write_representation(status)
       CALL channel_halt(substr, status)
    END IF

    CALL write_representation_dc(status, p_pe)
    CALL channel_halt(substr, status)

  END SUBROUTINE main_channel_initialize_reprs
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE associate_streams_to_channels

    ! ECHAM5
    USE mo_memory_base, ONLY: ostreams, nstreams, GAUSSIAN, SPECTRAL &
                            , FOURIER                                &
                            , t_stream, memory_info
    USE mo_linked_list, ONLY: list_element
    ! MESSy
    USE messy_main_grid_def_mem_bi, ONLY: nlev, nlevp1
    USE messy_main_mpi_bi,          ONLY: p_parallel_io
    USE messy_main_channel,         ONLY: new_channel, new_channel_object &
                                        , new_attribute                   &
                                        , AF_RST_CMP
    USE messy_main_channel_repr,    ONLY: t_representation, get_representation

    IMPLICIT NONE

    ! LOCAL
    INTEGER                      :: i
    TYPE (t_stream),     POINTER :: stream
    TYPE (list_element), POINTER :: element
    TYPE (list_element), TARGET  :: element_1st
    TYPE (memory_info),  POINTER :: info
    !
    CHARACTER(LEN=*), PARAMETER  :: substr = 'associate_streams_to_channels'
    INTEGER                      :: status
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: mem
    INTEGER                      :: reprid
    LOGICAL                      :: lok
    TYPE(t_representation), POINTER :: repr
    CHARACTER(len=32)            :: grid_type

    CALL start_message_bi(modstr,'STREAMS TO CHANNELS',substr)

    stream_loop: DO i=1, nstreams

       stream => ostreams(i)

       ! SKIP
       IF (TRIM(stream%name) == 'ls') CYCLE

       ! ---------------------------
       ! STREAM -> CHANNEL
       ! ---------------------------
       CALL new_channel(status, TRIM(stream%name))
       CALL channel_halt(substr, status)
       ! ---------------------------

       element_1st%next_list_element => stream%first_list_element
       element => element_1st
       element_loop: DO ! loop over elements in linked list
          element => element%next_list_element
          IF (.NOT.ASSOCIATED(element)) EXIT

          info => element%field%info
          grid_type = ''
          ! ---------------------------
          ! ELEMENT -> OBJECT
          ! ---------------------------
          SELECT CASE(info%repr)
          CASE(GAUSSIAN)
             !
             SELECT CASE(info%ndim)
                !
             CASE(1)
                lok = .FALSE.
             CASE(2)
                lok = .TRUE.
                reprid = GP_2D_HORIZONTAL
             CASE(3)
                IF (info%klev == 1) THEN
                   lok = .TRUE.
                   reprid = GP_3D_1LEV
                ELSE IF (info%klev == 2) THEN
                   lok = .TRUE.
                   reprid = GP_3D_2LEV
                ELSE IF (info%klev == 5) THEN
                   lok = .TRUE.
                   reprid = GP_3D_BELOWSF
                ELSE IF (info%klev == nlev) THEN
                   lok = .TRUE.
                   reprid = GP_3D_MID
                ELSE IF (info%klev == nlevp1) THEN
                   lok = .TRUE.
                   reprid = GP_3D_INT
                ELSE
                   lok = .FALSE.
                END IF
             CASE(4)
                lok = .FALSE.
             CASE DEFAULT
                lok = .FALSE.
                !
             END SELECT
             !
             grid_type = 'gaussian'
             !
          CASE(SPECTRAL)
             !
             SELECT CASE(info%ndim)
                !
             CASE(1)
                lok = .FALSE.
             CASE(2)
                IF ( (SIZE(element%field%ptr(:,:,:,:),2) == 2) .AND. & ! COMPLEX
                     (info%klev == 1) ) THEN
                   lok = .TRUE.
                   reprid = SP_2D_HORIZONTAL
                ELSE
                   lok = .FALSE.
                END IF
             CASE(3)
                IF (info%klev == 1) THEN
                   lok = .TRUE.
                   reprid = SP_2D_HORIZONTAL
                ELSE IF (info%klev == nlev) THEN
                   lok = .TRUE.
                   reprid = SP_3D_MID
                ELSE IF (info%klev == nlevp1) THEN
                   lok = .TRUE.
                   reprid = SP_3D_INT
                ELSE
                   lok = .FALSE.
                END IF
             CASE(4)
                lok = .FALSE.
             CASE DEFAULT
                lok = .FALSE.
                !
             END SELECT
             !
             grid_type = 'spectral, triangular truncation'
             !
          CASE(FOURIER)
             !
             SELECT CASE(info%ndim)
                !
             CASE(1)
                lok = .FALSE.
             CASE(2)
                lok = .TRUE.
                reprid = FAS_MID_ZM
             CASE(3)
                lok = .FALSE.
             CASE(4)
                IF (info%klev == 1) THEN
                   lok = .TRUE.
                   reprid = FAS_MID_ZM
                ELSE IF (info%klev == nlev) THEN
                   lok = .TRUE.
                   reprid = FAS_MID
                ELSE IF (info%klev == nlevp1) THEN
                   lok = .TRUE.
                   reprid = FAS_INT
                ELSE
                   lok = .FALSE.
                END IF
             CASE DEFAULT
                lok = .FALSE.
                !
             END SELECT
             !
          CASE DEFAULT
             !
             lok = .FALSE.
             !
          END SELECT

          IF (lok) THEN

             ! SPECIAL HANDLING FOR ECHAM5-STREAM-ELEMENTS WITH laccu = .TRUE.
             ! -> DOWNGRADE TO INSTANTANEOUS VALUES
             IF (info%laccu) THEN
                IF (p_parallel_io) THEN
                   WRITE(*,*) 'WARNING: stream/element ''' &
                        ,TRIM(stream%name),'/', TRIM(info%name) &
                        ,''' has laccu = .TRUE. !'
                END IF
                naccu = naccu + 1
             END IF

             mem => element%field%ptr(:,:,:,:)

             CALL new_channel_object(status, TRIM(stream%name) &
                  , TRIM(info%name)                            &
                  , mem = mem                                  &
                  , reprid = reprid                            &
                  , lrestreq = .NOT. info%contnorest .OR. info%lrerun )
             CALL channel_halt(substr, status)

             IF (TRIM(info%longname) /= '') THEN
                CALL new_attribute(status &
                     , TRIM(stream%name), TRIM(info%name), 'long_name' &
                     , c=TRIM(info%longname), iflag=AF_RST_CMP)
                CALL channel_halt(substr, status)
             END IF
             IF (TRIM(info%units) /= '') THEN
                CALL new_attribute(status &
                     , TRIM(stream%name), TRIM(info%name), 'units' &
                     , c=TRIM(info%units), iflag=AF_RST_CMP)
                CALL channel_halt(substr, status)
             END IF

             CALL get_representation(status, reprid, repr)
             CALL channel_halt(substr, status)
             CALL new_attribute(status &
                  , TRIM(stream%name), TRIM(info%name), 'representation' &
                  , c=TRIM(repr%name), iflag=AF_RST_CMP)
             CALL channel_halt(substr, status)

             ! additional attributes
             IF (TRIM(grid_type) /= '') THEN
                CALL new_attribute(status &
                     , TRIM(stream%name), TRIM(info%name), 'grid_type' &
                     , c=TRIM(grid_type))
                CALL channel_halt(substr, status)
             END IF

             IF (info%gribtable /= 0) THEN
                CALL new_attribute(status &
                     , TRIM(stream%name), TRIM(info%name), 'table' &
                     , i=info%gribtable)
                CALL channel_halt(substr, status)
             ENDIF

             IF (info%gribcode /= 0) THEN
                CALL new_attribute(status &
                     , TRIM(stream%name), TRIM(info%name), 'code' &
                     , i=info%gribcode)
                CALL channel_halt(substr, status)
             ENDIF

          ELSE

             IF (p_parallel_io) THEN
                WRITE(*,*) 'WARNING: stream/element ''' &
                     ,TRIM(stream%name),'/', TRIM(info%name) &
                     ,''' cannot be associated to channel!'
             END IF

          END IF
          ! ---------------------------

       END DO element_loop

    END DO stream_loop

    CALL end_message_bi(modstr,'STREAMS TO CHANNELS',substr)

  END SUBROUTINE associate_streams_to_channels
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE reset_accu_stream_elements(flag)

    ! ECHAM5
    USE mo_memory_base,     ONLY: ostreams, nstreams, t_stream, memory_info
    USE mo_linked_list,     ONLY: list_element
    USE messy_main_timer,   ONLY: delta_time

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)          :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'reset_accu_stream_elements'
    LOGICAL, SAVE                :: lfirst = .TRUE.
    INTEGER                      :: i, count
    TYPE (t_stream),     POINTER :: stream
    TYPE (list_element), POINTER :: element
    TYPE (list_element), TARGET  :: element_1st
    TYPE (memory_info),  POINTER :: info

    IF (lfirst) THEN
       ALLOCATE(ptr_accu(naccu))
       count = 0

       stream_loop: DO i=1, nstreams

          stream => ostreams(i)

          element_1st%next_list_element => stream%first_list_element
          element => element_1st
          element_loop: DO ! loop over elements in linked list
             element => element%next_list_element
             IF (.NOT.ASSOCIATED(element)) EXIT

             info => element%field%info

             ! ---------------------------
             IF (info%laccu) THEN
                count = count + 1
                ptr_accu(count)%ptr => element%field%ptr(:,:,:,:)
             END IF
             ! ---------------------------

          END DO element_loop

       END DO stream_loop

       lfirst = .FALSE.
    END IF

    SELECT CASE(flag)
    CASE(1)
       ! BEFORE OUTPUT ...
       ! DOWNGRADE ALL laccu-FIELDS TO INSTANTANEOUS VALUES
       ! -> REPLACE a = 0 + b*dt BY a = 0 + b
       DO count=1, naccu
          ptr_accu(count)%ptr(:,:,:,:) = &
               ptr_accu(count)%ptr(:,:,:,:) / delta_time
       END DO
    CASE(2)
       ! AFTER OUTPUT ...
       ! DOWNGRADE ALL laccu-FIELDS TO INSTANTANEOUS VALUES
       ! -> REPLACE a = a + b*dt BY a = 0 + b*dt
       DO count=1, naccu
          ptr_accu(count)%ptr(:,:,:,:) = 0.0_DP
       END DO
    END SELECT

  END SUBROUTINE reset_accu_stream_elements
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE timer_sync

    USE mo_time_control,       ONLY: resume_date, start_date
    USE mo_time_conversion,    ONLY: TC_set, time_native, TC_convert
    USE messy_main_timer_bi,   ONLY: timer_message
    USE messy_main_timer,      ONLY: timer_get_date

    IMPLICIT NONE

    ! synchronize ECHAM5 timer with MESSy timer

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'timer_sync'
    TYPE(time_native)           :: io_date
    INTEGER                     :: yr, mo, dy, hr, mi, se, status

    CALL timer_get_date(status, "start", yr, mo, dy, hr, mi, se)
    CALL timer_message(status, substr)
    CALL TC_set(yr, mo, dy, hr, mi, se, io_date)
    CALL TC_convert(io_date, start_date)

    CALL timer_get_date(status, "resume", yr, mo, dy, hr, mi, se)
    CALL timer_message(status, substr)
    CALL TC_set(yr, mo, dy, hr, mi, se, io_date)
    CALL TC_convert(io_date, resume_date)

  END SUBROUTINE timer_sync
  ! -------------------------------------------------------------------
