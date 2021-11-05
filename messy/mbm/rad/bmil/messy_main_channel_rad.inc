! -*- f90 -*-
  ! -------------------------------------------------------------------
  !  CHANNEL INTERFACE FOR MBM OUTPUT
  ! -------------------------------------------------------------------

  !PRIVATE SUBROUTINE main_channel_initialize_gatts
  !PRIVATE SUBROUTINE main_channel_initialize_dims
  !PRIVATE SUBROUTINE main_channel_initialize_reprs
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_initialize_gatts

    USE messy_main_data_bi,       ONLY: basemodver=>modver, basemodstr=>modstr

    USE messy_main_channel,       ONLY: new_attribute, write_attribute &
                                      , AF_RST_CMP, AF_RST_INP, EXP_NAME &
                                      , modstr, modver, EXEC_CHECKSUM
    ! op_pj_20130314+
    USE messy_main_compilerinfo_mem, ONLY: compiler_version, compiler_call &
                                         , compiler_flags &
                                         , compiler_cppdefs, compiler_includes
    USE messy_main_constants_mem,    ONLY: STRLEN_ULONG
    ! op_pj_20130314-

    IMPLICIT NONE

    INTRINSIC :: DATE_AND_TIME, TRIM, LEN_TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_initialize_gatts'
    INTEGER :: status
    CHARACTER (8)   :: ydate
    CHARACTER (10)  :: ytime
    ! op_pj_20110222+
    CHARACTER (LEN=2)            :: istr
    INTEGER :: i1, i2, il, ic
    ! op_pj_20110222-

    ! example 1: system date and time
    CALL DATE_AND_TIME(ydate, ytime)
    CALL new_attribute(status, 'operating_date_time' &
         , c = ydate(1:8)//' '//ytime(1:6) )
    CALL channel_halt(substr, status)

    ! example 2: experiment name from namelist
    CALL new_attribute(status, 'MESSy_experiment', c=TRIM(EXP_NAME))
    CALL channel_halt(substr, status)

    ! CHECKSUM OF EXECUTABLE
    CALL new_attribute(status, 'EXEC_CHECKSUM', c=TRIM(EXEC_CHECKSUM))
    CALL channel_halt(substr, status)

    ! example 3: basemodel version information
    CALL new_attribute(status, 'MESSy_basemodel' &
         , c = basemodstr//' version '//basemodver//&
         &', IPA Uni Mainz and DLR-IPA Oberpfaffenhofen' )
    CALL channel_halt(substr, status)

    ! example 4: special flag to be checked from restart file
    CALL new_attribute(status, 'restart_test_flag' &
         , i=31415, iflag=AF_RST_CMP)
    CALL channel_halt(substr, status)

    ! op_pj_20130314+
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
    ! op_pj_20130314-

    ! diagnostic output 
    CALL write_attribute(status)
    CALL channel_halt(substr, status)

  END SUBROUTINE main_channel_initialize_gatts
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_initialize_dims

    USE messy_main_grid_def_mem_bi, ONLY: nvclev, vct &
                                        , nlev, nlon, nlat &
                                        , nglon, nglat     &
                                        , nmp1, lv_echam
    USE messy_main_grid_def_bi,     ONLY: philon, philat, pf, ph
    USE messy_main_timer,     ONLY: delta_time, current_time_step &
                                  , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND &
                                  , YEAR_START, MONTH_START, DAY_START     &
                                  , HOUR_START, MINUTE_START, SECOND_START

    USE messy_main_channel_dimensions,   ONLY: new_dimension            &
                                             , write_dimension          &
                                             , add_dimension_variable   &
                                             , add_dimension_variable_att

    IMPLICIT NONE
    INTRINSIC :: RANDOM_SEED, ABS, REAL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_initialize_dims'
    INTEGER                             :: status
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: array
    CHARACTER(LEN=80)                   :: start_date_string
    REAL(dp)                            :: yyyymmdd
    INTEGER                             :: n, i

    ! #############
    ! ### TIME ####
    ! #############
    ! ... dimension ...
    CALL new_dimension(status, DIMID_TIME, 'time', 1, .TRUE.)
    CALL channel_halt(substr, status)

    ! ... dimension variable 'time' ...
    CALL add_dimension_variable(status, 'time', 'time', (/ 0.0_dp /))
    CALL channel_halt(substr, status)
    ! ... with attribute 'long_name' ...
    CALL add_dimension_variable_att(status, 'time', 'time', &
         'long_name', c='time')
    CALL channel_halt(substr, status)
    ! ... with attribute 'units' ...

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

    ! ###############
    ! ### LEVELS ####
    ! ###############
    IF (lv_echam) THEN
       ! MID-LAYERS (FULL LEVELS)
       ! ------------------------
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
       !
       ! INTERFACE-LAYERS (HALF LEVELS)
       ! ------------------------------
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
    ELSE
       ! ########################
       ! ### PRESSURE LEVELS ####
       ! ########################

       ! MID-LAYERS (FULL PRESSURE LEVELS)
       !
       CALL new_dimension(status, DIMID_LEV, 'lev', nlev)
       CALL channel_halt(substr, status)
       !
       CALL add_dimension_variable(status, 'lev', 'lev', pf)
       CALL channel_halt(substr, status)
       !
       CALL add_dimension_variable_att(status, 'lev', 'lev', &
               'long_name', c='pressure level at layer midpoints')
       CALL channel_halt(substr, status)
       !
       CALL add_dimension_variable_att(status, 'lev', 'lev', &
               'standard_name', c='pressure')
       CALL channel_halt(substr, status)
       !
       CALL add_dimension_variable_att(status, 'lev', 'lev', &
               'units', c='Pa')
       CALL channel_halt(substr, status)
       !
       CALL add_dimension_variable_att(status, 'lev', 'lev', &
               'positive', c='down')
       CALL channel_halt(substr, status)
       !
       CALL add_dimension_variable_att(status, 'lev', 'lev', &
               'formula', c='hyam hybm (press=hyam+hybm*aps)')
       CALL channel_halt(substr, status)

       ! INTERFACE-LAYERS (HALF PRESSURE LEVELS)
       !
       CALL new_dimension(status, DIMID_ILEV, 'ilev', nlev+1)
       CALL channel_halt(substr, status)
       !
       CALL add_dimension_variable(status, 'ilev', 'ilev', ph)
       CALL channel_halt(substr, status)
       !
       CALL add_dimension_variable_att(status, 'ilev', 'ilev', &
               'long_name', c='pressure level at layer interfaces')
       CALL channel_halt(substr, status)
       !
       CALL add_dimension_variable_att(status, 'ilev', 'ilev', &
               'standard_name', c='pressure')
       CALL channel_halt(substr, status)
       !
       CALL add_dimension_variable_att(status, 'ilev', 'ilev', &
               'units', c='Pa')
       CALL channel_halt(substr, status)
       !
       CALL add_dimension_variable_att(status, 'ilev', 'ilev', &
               'positive', c='down')
       CALL channel_halt(substr, status)
       !
       CALL add_dimension_variable_att(status, 'ilev', 'ilev', &
               'formula', c='hyai hybi (ilev=hyai+hybi*aps)')
       CALL channel_halt(substr, status)

    END IF
    
    ! ###################
    ! ### LONGITUDES ####
    ! ###################

    ! ... dimension ...
    CALL new_dimension(status, DIMID_LON, 'lon', nlon)
    CALL channel_halt(substr, status)

    ! ... dimension variable ...
    ALLOCATE(array(nlon))
    array = philon
    CALL add_dimension_variable(status, 'lon', 'lon', array)
    CALL channel_halt(substr, status)
    DEALLOCATE(array)

    ! ... with attribute 'long_name' ...
    CALL add_dimension_variable_att(status, 'lon', 'lon', &
         'long_name', c='longitude')
    CALL channel_halt(substr, status)

    ! ... with attribute 'units' ...
    CALL add_dimension_variable_att(status, 'lon', 'lon', &
         'units', c='degrees_east')
    CALL channel_halt(substr, status)

    ! ##################
    ! ### LATITUDES ####
    ! ##################

    ! ... dimension ...
    CALL new_dimension(status, DIMID_LAT, 'lat', nlat)
    CALL channel_halt(substr, status)

    ! ... dimension variable ...
    ALLOCATE(array(nlat))
    array = philat
    CALL add_dimension_variable(status, 'lat', 'lat', array)
    CALL channel_halt(substr, status)
    DEALLOCATE(array)

    ! ... with attribute 'long_name' ...
    CALL add_dimension_variable_att(status, 'lat', 'lat', &
         'long_name', c='latitude')
    CALL channel_halt(substr, status)

    ! ... with attribute 'units' ...
    CALL add_dimension_variable_att(status, 'lat', 'lat', &
         'units', c='degrees_north')
    CALL channel_halt(substr, status)

    ! ###############################
    ! ### ARBITRARY ARRAY LENGTH ####
    ! ###############################
    
    ! Note: Here the length of the random state vector of the fortran system
    !       is arbitrarily used as array length ...
    CALL RANDOM_SEED(SIZE = n)

    ! ... dimension ...
    CALL new_dimension(status, DIMID_AL, 'n', n)
    CALL channel_halt(substr, status)
    
    ! ################
    ! ### SPECIAL ####
    ! ################

    ! 2 boundaries
    CALL new_dimension(status, DIMID_TBNDS, 'tbnds', tbnds)
    CALL channel_halt(substr, status)

    ! #################
    ! DIAGNOSTIC OUTPUT
    ! #################
    CALL write_dimension(status)
    CALL channel_halt(substr, status)

  END SUBROUTINE main_channel_initialize_dims
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_initialize_reprs

    USE messy_main_channel_repr,         ONLY: new_representation        &
                                             , write_representation_dc   &
                                             , set_representation_decomp &
                                             , write_representation      &
                                             , AUTO, IRANK               &
                                             , PIOTYPE_SGL, PIOTYPE_IND  &
                                             , PIOTYPE_COL
    USE messy_main_channel_dimensions,   ONLY: new_dimension
    USE messy_main_grid_def_mem_bi,      ONLY: glons, glats  &
                                             , glone, glate  &
                                             , nproma, ngpblks, nlev &
                                             , nglon, nglat

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'main_channel_initialize_reprs'
    INTEGER                        :: status
    INTEGER, DIMENSION(:), POINTER :: scdim => NULL()
    !
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg    = 0  ! max no. of segments
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()  ! start vector
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()  ! count vector
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()  ! memory lower bound
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()  ! memory upper bound
    INTEGER                          :: i
    
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

    CALL new_representation(status, GP_3D_MID, 'GP_3D_MID'        &
         , rank = 3, link = 'xxx-', dctype = DC_GP                &
         , dimension_ids = (/ DIMID_LON, DIMID_LAT, DIMID_LEV /)  &
         , ldimlen       = (/ nproma   , ngpblks,   AUTO      /)  &
         , output_order  = (/ 1,2,3 /)                            &
         , axis = 'XYZ-'                                          &
         )
    CALL channel_halt(substr, status)

    start(:,1) = glons(:)
    start(:,2) = glats(:)
    start(:,3) = 1
    
    cnt(:,1) = glone(:)-glons(:)+1
    cnt(:,2) = glate(:)-glats(:)+1
    cnt(:,3) = nlev

    meml(:,1) = 1
    !meml(:,2) = (/ 1, nglh(1)+1 /)
    meml(:,2) = 1
    meml(:,3) = 1

    memu(:,1) = nglon
    !memu(:,2) = (/ nglh(1), nglh(1)+nglh(2)/)
    memu(:,2) = nglat
    memu(:,3) = nlev
    
    DO i = 1, IRANK
       PRINT *,'start(',i,') = ',start(:,i)
       PRINT *,'cnt  (',i,') = ',cnt  (:,i)
       PRINT *,'meml (',i,') = ',meml (:,i)
       PRINT *,'memu (',i,') = ',memu (:,i)
       PRINT *
    END DO
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

    gp_start(:,3) = start(:,3)
    gp_cnt(:,3) = cnt(:,3)
    gp_meml(:,3) = meml(:,3)
    gp_memu(:,3) = memu(:,3)

    gp_start(:,4) = 1
    gp_cnt(:,4) = 1
    gp_meml(:,4) = 1
    gp_memu(:,4) = 1

    ! -----------------------------------------------------

    CALL new_representation(status, GP_3D_INT, 'GP_3D_INT'         &
         , rank = 3, link = 'xxx-', dctype = DC_GP                 &
         , dimension_ids = (/ DIMID_LON, DIMID_LAT, DIMID_ILEV /)  &
         , ldimlen       = (/ nproma   , ngpblks,   AUTO   /) &
         , output_order  = (/ 1,2,3 /)                             &
         , axis = 'XYZ-'                                           &
         )
    CALL channel_halt(substr, status)

    cnt(:,3)  = nlev+1
    memu(:,3) = nlev+1

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

    !start(:,2) = start(:,2)
    !cnt(:,2)   = cnt(:,3)
    !meml(:,2)  = meml(:,3)
    !memu(:,2)  = memu(:,3)

    start(:,3) = 1
    cnt(:,3)   = 1
    meml(:,3)  = 1
    memu(:,3)  = 1

    CALL set_representation_decomp(status, GP_2D_HORIZONTAL &
        , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -----------------------------------------------------

    CALL new_representation(status, GP_3D_1LEV, 'GP_3D_1LEV' &
         , rank = 2, link = 'xx--', dctype = DC_GP              &
         , dimension_ids = (/ DIMID_LON, DIMID_LAT /)    &
         , ldimlen       = (/ nproma   , ngpblks   /)    &
         , axis = 'XY--'                                 &
         )
    CALL channel_halt(substr, status)
    
    ! NOTE: THE RANK IS W.R.T. MEMORY, I.E., ACCORDING TO 'link' !!!

    !start(:,3) = start(:,2)
    !cnt(:,3)   = cnt(:,2)
    !meml(:,3)  = meml(:,2)
    !memu(:,3)  = memu(:,2)

    !start(:,2) = 1
    !cnt(:,2)   = 1
    !meml(:,2)  = 1
    !memu(:,2)  = 1

    CALL set_representation_decomp(status, GP_3D_1LEV &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
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

    ! ################################
    ! ### GRIDPOINT WITH BOUNDARY ####
    ! ################################
    !!$CALL new_representation(status, GP_3D_MID_BND, 'GP_3D_MID_BND'   &
    !!$     , rank = 3, link = 'xxx-', dctype = 0                       &
    !!$     , dimension_ids = (/ DIMID_LON, DIMID_LAT, DIMID_LEV /)     &
    !!$     , ldimlen       = (/ AUTO     , AUTO     , AUTO      /)     &
    !!$     , nbounds       = (/ 2,2,0 /)                               &
    !!$     , output_order  = (/ 1,2,3 /)                               &
    !!$     , axis = 'XYZ-'                                             &
    !!$     )
    !!$CALL channel_halt(substr, status) 
    
    ! ###############
    ! ### SCALAR ####
    ! ###############
    ALLOCATE(scdim(0))
    CALL new_representation(status, SCALAR, 'SCALAR'                 &
         , rank = 0, link = '----', dctype = 0                       &
         , dimension_ids = scdim &
         , ldimlen       = scdim &
         , axis = '----'                                             &
         )
    CALL channel_halt(substr, status)
    DEALLOCATE(scdim)
    NULLIFY(scdim)

    ! ##############
    ! ### ARRAY ####
    ! ##############
    !!$CALL new_representation(status, ARRAY, 'ARRAY'                 &
    !!$     , rank = 1, link = 'x---', dctype = 0                     &
    !!$     , dimension_ids = (/ DIMID_AL /) &
    !!$     , ldimlen       = (/ AUTO /)     &
    !!$     , axis = 'N---'                                           &
    !!$     )
    !!$CALL channel_halt(substr, status)

    ! #################
    ! DIAGNOSTIC OUTPUT
    ! #################
    CALL write_representation(status)
    CALL channel_halt(substr, status)

  END SUBROUTINE main_channel_initialize_reprs
  ! -------------------------------------------------------------------
  ! -------------------------------------------------------------------
