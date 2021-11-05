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
         &', MPIC, Mainz' )
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
    CHARACTER(LEN=80)                   :: start_date_string
    REAL(dp)                            :: yyyymmdd

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
    CALL add_dimension_variable_att(status, 'time', 'time', &
         'bounds', c='time_bnds')
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

    ! ################
    ! ### SPECIAL ####
    ! ################
    ! time boundaries
    CALL new_dimension(status, DIMID_TBNDS, 'tbnds', 2)
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

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'main_channel_initialize_reprs'
    INTEGER                        :: status
    INTEGER, DIMENSION(:), POINTER :: scdim => NULL()

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

    ! ####################
    ! ### TIME BOUNDS ####
    ! ####################
    CALL new_representation(status, REPR_TIMEBNDS, 'TIME_BNDS' &
         , rank = 1, link = 'x---', dctype = DC_BC             &
         , dimension_ids = (/DIMID_TBNDS/)                     &
         , ldimlen       = (/ tbnds     /)                     &
         , axis = 'N---'                                       &
         )
    CALL channel_halt(substr, status)

    ! #################
    ! DIAGNOSTIC OUTPUT
    ! #################
    CALL write_representation(status)
    CALL channel_halt(substr, status)

  END SUBROUTINE main_channel_initialize_reprs
  ! -------------------------------------------------------------------
