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
#ifdef MBM_CLAMS
  SUBROUTINE main_channel_initialize_gatts

#else
  SUBROUTINE clams_initialize_gatts
#endif
    use netcdf
    
    USE messy_main_data_bi,       ONLY: basemodver=>modver, basemodstr=>modstr

    USE messy_main_channel,       ONLY: new_attribute, write_attribute &
                                      , AF_RST_CMP, AF_RST_INP, EXP_NAME &
                                      , modstr, modver
    USE messy_main_mpi_bi,        ONLY: p_pe

    USE messy_clams_global,       ONLY: ldiagout, username, initfile

#ifndef MBM_CLAMS
    USE messy_main_channel_error_bi,  ONLY: channel_halt
    USE messy_main_constants_mem,     ONLY: DP
#endif

    IMPLICIT NONE

    INTRINSIC :: DATE_AND_TIME, TRIM

    ! LOCAL
#ifdef MBM_CLAMS
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_initialize_gatts'
#else
    CHARACTER(LEN=*), PARAMETER :: substr = 'clams_initialize_gatts'
#endif
    CHARACTER(nf90_max_name) :: attname, cval
    INTEGER  :: ncid, natts, iatt, xtype, ival
    REAL(DP) :: rval

    INTEGER :: status
    CHARACTER (8)   :: ydate
    CHARACTER (10)  :: ytime
    CHARACTER (30)  :: timestr

!#ifdef MBM_CLAMS
    CALL new_attribute(status, 'pi_name', c = username)
    CALL channel_halt (substr, status)
    CALL new_attribute(status, 'organisation', c='Forschungszentrum Juelich GmbH')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'institut', c='IEK-7')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'address', c='52425 Juelich')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'data_category', c='MODEL')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'experiment', c='CLaMS')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'platform', c='Computer')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'type_of_data', c='SIMULATION')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'data_origin', c='ECMWF')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'data_description', c='GRIB')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'title', c='3D Trajectories')
    CALL channel_halt(substr, status)
    CALL DATE_AND_TIME(ydate, ytime)
    WRITE(timestr,*)   &
         ydate(1:4),'-',ydate(5:6),'-',ydate(7:8),' ',  &
         ytime(1:2),':',ytime(3:4),':',ytime(5:6)
    CALL new_attribute(status, 'file_creation_time' &
         , c = timestr )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'file_modification_time' &
         , c = timestr )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'history', &
         c='Created on '//trim(timestr)//' MESSy-CLaMS')
    CALL channel_halt(substr, status)

    ! copy attributes  exp_POS_xxx and exp_VERTCOOR_xxx from initfile
    status = nf90_open (initfile, nf90_nowrite, ncid)
    if (status==0) then
       status = nf90_inquire (ncid, nAttributes=natts)
       do iatt = 1, natts
          status = nf90_inq_attname (ncid, nf90_global, iatt, attname)
          if (attname(1:7)=='exp_POS' .or. attname(1:12)=='exp_VERTCOOR') then
             status = nf90_inquire_attribute (ncid, nf90_global, attname, xtype=xtype)
             SELECT CASE (xtype)
             CASE (NF90_INT)
                status = nf90_get_att (ncid, nf90_global, attname, ival)
                CALL new_attribute(status, trim(attname), i=ival)
                CALL channel_halt(substr, status)                
             CASE (NF90_FLOAT, NF90_DOUBLE)
                status = nf90_get_att (ncid, nf90_global, attname, rval)
                CALL new_attribute(status, trim(attname), r=rval)
                CALL channel_halt(substr, status)                
             CASE (NF90_CHAR)
                status = nf90_get_att (ncid, nf90_global, attname, cval)
                CALL new_attribute(status, trim(attname), c=trim(cval))
                CALL channel_halt(substr, status)          
             END SELECT
           endif
       enddo
    endif
    status = nf90_close (ncid)


!#endif

    ! example 3: basemodel version information
    CALL new_attribute(status, 'MESSy_basemodel' &
         , c = basemodstr//' version '//basemodver//&
         &', IPA Uni Mainz and DLR-IPA Oberpfaffenhofen' )
    CALL channel_halt(substr, status)

    ! example 4: special flag to be checked from restart file
    CALL new_attribute(status, 'restart_test_flag' &
         , i=31415, iflag=AF_RST_CMP)
    CALL channel_halt(substr, status)

    ! diagnostic output 
    if (p_pe==0 .and. ldiagout) then
       CALL write_attribute(status)
       CALL channel_halt(substr, status)
    endif

#ifdef MBM_CLAMS
  END SUBROUTINE main_channel_initialize_gatts
#else
  END SUBROUTINE clams_initialize_gatts
#endif
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
#ifdef MBM_CLAMS
  SUBROUTINE main_channel_initialize_dims
#else
  SUBROUTINE clams_initialize_dims
#endif

    USE messy_main_channel_dimensions,   ONLY: new_dimension            &
                                             , write_dimension          &
                                             , add_dimension_variable   &
                                             , add_dimension_variable_att
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check,  &
                                        read_nml_close, find_next_free_unit
    USE messy_clams_global,       ONLY: ntasks, nparts_max, dnparts_max, & 
                                        nparts_max_shuffle,&
                                        nx, ny, nz, ntheta, &
                                        longrid, latgrid, levelgrid, thetagrid
    USE messy_clamsmix_global,    ONLY: adapt_par
    USE messy_cmn_photol_mem,     ONLY: IP_MAX
!!$    USE messy_dissoc,             ONLY: jpslev, jpslevall, jpwave, jpschi,&
!!$                                        jplats, numj_max
    USE messy_main_timer,         ONLY: delta_time, current_time_step, &
                                        YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, &
                                        YEAR_START, MONTH_START, DAY_START,     &
                                        HOUR_START, MINUTE_START, SECOND_START
    USE messy_clamssedi_global,   ONLY: nparticle_max

!!$#ifndef MBM_CLAMS
!!$    USE messy_main_channel_bi,    ONLY: DIMID_TRAJ, DIMID_NTASKS,&
!!$         DIMID_SHUFFLE, channel_halt, DIMID_WAVE, DIMID_DISSOC_LEV, &
!!$         DIMID_DISSOC_LEVC, DIMID_DISSOC_LAT, DIMID_DISSOC_SZA, DIMID_DISSOC_NUMJ, &
!!$         DIMID_SEDI_PARTICLE, DIMID_MIX_GRID
!!$    USE messy_main_constants_mem, ONLY: DP
!!$#endif
#ifndef MBM_CLAMS
    USE messy_main_channel_error_bi, ONLY: channel_halt 
    USE messy_main_channel_bi,       ONLY: DIMID_TRAJ, DIMID_NTASKS,&
                                           DIMID_SHUFFLE,  &
                                           DIMID_SEDI_PARTICLE, DIMID_MIX_GRID, &
                                           DIMID_TIME
    USE messy_main_constants_mem,    ONLY: DP
#endif

    IMPLICIT NONE
   
    ! LOCAL
#ifdef MBM_CLAMS
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_initialize_dims'
#else
    CHARACTER(LEN=*), PARAMETER :: substr = 'clams_initialize_dims'
#endif
    INTEGER                     :: status

    CHARACTER(LEN=80)           :: start_date_string
    REAL(DP)                    :: yyyymmdd
    INTEGER                     :: i

    ! #############
    ! ### TIME ####
    ! #############
    ! ... dimension ...
#ifdef MBM_CLAMS
    CALL new_dimension(status, DIMID_TIME, 'time', 1, .TRUE.)
    CALL channel_halt(substr, status)
#endif

    ! Dimension variable JULIAN_SECONDS
    CALL add_dimension_variable(status, 'time', 'JULIAN_SECONDS', (/ 0.0_dp /))
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'time', 'JULIAN_SECONDS', &
         'long_name', c='Time')
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'time', 'JULIAN_SECONDS', &
         'units', c='seconds since 2000-01-01 00:00:00 UTC')
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'time', 'JULIAN_SECONDS', &
         'description', c='Time in julian seconds')
    CALL channel_halt(substr, status)
    
#ifdef MBM_CLAMS
    ! ... dimension variable 'time' ...
    CALL add_dimension_variable(status, 'time', 'time', (/ 0.0_dp /))
    CALL channel_halt(substr, status)
    ! ... with attribute 'long_name' ...
    CALL add_dimension_variable_att(status, 'time', 'time', &
         'long_name', c='Time')
    ! ... with attribute 'time_bnds' ...
    CALL add_dimension_variable_att(status, 'time', 'time', &
         'bounds', c='time_bnds')
    CALL channel_halt(substr, status)
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

    ! ##################
    ! ###  PARCELS  ####
    ! ##################
    
    CALL new_dimension(status, DIMID_LON, 'LON', nx)
    CALL channel_halt(substr, status)
    CALL add_dimension_variable(status,'LON','LON',longrid)
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'LON', 'LON', 'units', c='deg E')
    CALL channel_halt(substr, status)

    CALL new_dimension(status, DIMID_LAT, 'LAT', ny)
    CALL channel_halt(substr, status)
    CALL add_dimension_variable(status,'LAT','LAT',latgrid)
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'LAT', 'LAT', 'units', c='deg N')
    CALL channel_halt(substr, status)

    CALL new_dimension(status, DIMID_LEV, 'LEV', nz)
    CALL channel_halt(substr, status)
    CALL add_dimension_variable(status,'LEV','LEV',levelgrid(1:nz))
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'LEV', 'LEV', 'units', c=' ')
    CALL channel_halt(substr, status)

    if (ntheta>0) then
       CALL new_dimension(status, DIMID_THETA, 'THETA', ntheta)
       CALL channel_halt(substr, status)
       CALL add_dimension_variable(status,'THETA','THETA',thetagrid(1:ntheta))
       CALL channel_halt(substr, status)
       CALL add_dimension_variable_att(status, 'THETA', 'THETA', 'units', c='K')
       CALL channel_halt(substr, status)
    endif
#endif

    CALL new_dimension(status, DIMID_TRAJ, 'NPARTS', nparts_max)
    CALL channel_halt(substr, status)

    CALL new_dimension(status, DIMID_SHUFFLE, 'NPARTS_SHUFFLE', nparts_max_shuffle)
    CALL channel_halt(substr, status)

    CALL new_dimension(status, DIMID_MIX_GRID, 'NLEVS', adapt_par%nlevs)
    CALL channel_halt(substr, status)

!!$! jug define DISSOC dimensions
!!$    CALL new_dimension(status, DIMID_WAVE, 'WAVELENGTH', jpwave)
!!$    CALL channel_halt(substr, status)
!!$
!!$    CALL new_dimension(status, DIMID_DISSOC_LEV, 'DISSOC_LEV', jpslevall)
!!$    CALL channel_halt(substr, status)
!!$
!!$    CALL new_dimension(status, DIMID_DISSOC_LEVC, 'DISSOC_LEVC', jpslev)
!!$    CALL channel_halt(substr, status)
!!$
!!$    CALL new_dimension(status, DIMID_DISSOC_LAT, 'DISSOC_LAT', jpschi)
!!$    CALL channel_halt(substr, status)
!!$
!!$    CALL new_dimension(status, DIMID_DISSOC_SZA, 'DISSOC_SZA', jpschi)
!!$    CALL channel_halt(substr, status)
!!$
!!$    ! set number of photolysis rates -- to be replaced by smaller number 
!!$    ! if determined before
!!$    CALL new_dimension(status, DIMID_DISSOC_NUMJ, 'DISSOC_NUMJ', numj_max)
!!$    CALL channel_halt(substr, status)
!!$
    ! define SEDI dimensions
    CALL new_dimension(status, DIMID_SEDI_PARTICLE, 'SEDI_PARTICLE', nparticle_max)
    CALL channel_halt(substr, status)


!!!!!
    CALL new_dimension(status, DIMID_NTASKS, 'ntasks', ntasks)
    CALL channel_halt(substr, status)

#ifdef MBM_CLAMS
    ! ################
    ! ### SPECIAL ####
    ! ################

    ! 2 boundaries
     CALL new_dimension(status, DIMID_BNDS, 'bnds', 2)
     CALL channel_halt(substr, status)
#endif
  

    
#ifdef MBM_CLAMS
  END SUBROUTINE main_channel_initialize_dims
#else
  END SUBROUTINE clams_initialize_dims
#endif
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
#ifdef MBM_CLAMS
  SUBROUTINE main_channel_initialize_reprs
#else
  SUBROUTINE clams_initialize_reprs
#endif

#ifndef MBM_CLAMS
    USE messy_main_channel_bi,           ONLY: DC_IX
#endif
    USE messy_main_mpi_bi,               ONLY: p_pe
    USE messy_clams_global,              ONLY: dnparts_max, dnparts_max_shuffle,&
                                               ldiagout, ntheta
    USE messy_clamsmix_global,           ONLY: adapt_par
!!$    USE messy_dissoc,                    ONLY: jpslev, jpslevall, jpwave,&
!!$                                               jpschi, jplats, numj_max
    USE messy_cmn_photol_mem,            ONLY: IP_MAX
    USE messy_main_channel_repr,         ONLY: new_representation &
         , set_representation_decomp, write_representation &
         , AUTO, IRANK               &
         , PIOTYPE_SGL, PIOTYPE_IND  &
         , PIOTYPE_COL
    USE messy_clamssedi_global,          ONLY: nparticle_max

!!$#ifndef MBM_CLAMS
!!$    USE messy_main_channel_bi, ONLY: REPR_LG_CLAMS, REPR_CLAMS_SHUFFLE, &
!!$         REPR_3DINP_CLAMS, REPR_3DINP_CLAMSTHETA, &
!!$         REPR_NTASKS,&
!!$         REPR_DISSOC_WAVE, REPR_DISSOC_LEV,  REPR_DISSOC_LEVC, &
!!$         REPR_DISSOC_SZA, REPR_DISSOC_LAT, REPR_DISSOC_TABS,  &
!!$         REPR_DISSOC_TABS_DAVG, REPR_DISSOC_TWOD, REPR_DISSOC_TWODC, &
!!$         REPR_DISSOC_RATE, REPR_SEDI_PARTICLE, REPR_MIX_GRID, &
!!$         DIMID_TRAJ, DIMID_NTASKS, DIMID_SHUFFLE, &
!!$         DIMID_LON, DIMID_LAT, DIMID_LEV, DIMID_THETA, channel_halt, &
!!$         DIMID_WAVE, DIMID_DISSOC_LEV, DIMID_DISSOC_LEVC, &
!!$         DIMID_DISSOC_LAT, DIMID_DISSOC_SZA, DIMID_DISSOC_NUMJ, &
!!$         DIMID_SEDI_PARTICLE, DIMID_MIX_GRID
!!$#endif
#ifndef MBM_CLAMS
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: REPR_LG_CLAMS, REPR_CLAMS_SHUFFLE, &
         REPR_3DINP_CLAMS, REPR_3DINP_CLAMSTHETA, &
         REPR_NTASKS,&
         REPR_SEDI_PARTICLE, REPR_MIX_GRID, &
         DIMID_TRAJ, DIMID_NTASKS, DIMID_SHUFFLE, &
         DIMID_LON, DIMID_LAT, DIMID_LEV, DIMID_THETA, &
         DIMID_SEDI_PARTICLE, DIMID_MIX_GRID
#endif

    IMPLICIT NONE 

    ! LOCAL
#ifdef MBM_CLAMS
    CHARACTER(LEN=*), PARAMETER    :: substr = 'main_channel_initialize_reprs'
#else
    CHARACTER(LEN=*), PARAMETER    :: substr = 'clams_initialize_reprs'
#endif
    INTEGER                        :: status
    INTEGER, DIMENSION(:), POINTER :: scdim => NULL()

    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg    = 0  ! max no. of segments
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    !------------------------------------------
    !  SCALAR 
    !------------------------------------------
#ifdef MBM_CLAMS
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
#endif

    !------------------------------------------
    ! For dnparts:
    !------------------------------------------
    CALL new_representation(status, REPR_NTASKS, 'REPR_NTASKS'       &
         , rank = 1, link = 'x---', dctype = DC_IX                   &
         , dimension_ids = (/ DIMID_NTASKS /) &
         , ldimlen       = (/ 1 /) &
         , axis = 'N---'                                             &
         )
    CALL channel_halt(substr, status)

    !------------------------------------------------
    ! For winddata (changed order of lon,lat,lev)
    !------------------------------------------------
    CALL new_representation(status, REPR_3DINP_CLAMS, 'REPR_3DINP_CLAMS' &
         , rank = 3, link = 'xxx-', dctype = 0                           &
         , dimension_ids = (/ DIMID_LON, DIMID_LAT, DIMID_LEV /)         &
         , ldimlen = (/ AUTO, AUTO, AUTO /)                              &
         , axis = 'XYZ-'                                                 &
         )
    CALL channel_halt(substr, status)

    !------------------------------------------------
    ! Parameters from met. files on theta levels !!!
    !------------------------------------------------
    if (ntheta>0) then
       CALL new_representation(status, REPR_3DINP_CLAMSTHETA, 'REPR_3DINP_CLAMSTHETA' &
            , rank = 3, link = 'xxx-', dctype = 0                           &
            , dimension_ids = (/ DIMID_LON, DIMID_LAT, DIMID_THETA /)         &
            , ldimlen = (/ AUTO, AUTO, AUTO /)                              &
            , axis = 'XYZ-'                                                 &
            )
       CALL channel_halt(substr, status)
    endif

    !----------------------------------------------------------
    ! clams parcels
    !----------------------------------------------------------
    CALL new_representation(status, REPR_LG_CLAMS, 'REPR_LG_CLAMS'         &
         , rank = 1, link = 'x---', dctype = DC_IX                 &
         , dimension_ids = (/ DIMID_TRAJ /) &
         , ldimlen       = (/ dnparts_max /) &
         , axis = 'N---'                                           &
         )
    CALL channel_halt(substr, status)

    CALL new_representation(status, REPR_CLAMS_SHUFFLE, 'REPR_CLAMS_SHUFFLE'  &
         , rank = 1, link = 'x---', dctype = DC_IX                       &
         , dimension_ids = (/ DIMID_SHUFFLE /)                              &
         , ldimlen       = (/ dnparts_max_shuffle /)                              &
         , axis = 'N---'                                                 &
         )
    CALL channel_halt(substr, status)

    CALL new_representation(status, REPR_MIX_GRID, 'REPR_MIX_GRID'         &
         , rank = 1, link = 'x---', dctype = 0                     &
         , dimension_ids = (/ DIMID_MIX_GRID /) &
         , ldimlen       = (/ adapt_par%nlevs /) &
         , axis = 'N---'                                           &
         )
    CALL channel_halt(substr, status)

    
    !----------------------------------------------------------
    ! SEDI representations
    !----------------------------------------------------------
    CALL new_representation(status, REPR_SEDI_PARTICLE, 'REPR_SEDI_PARTICLE'  &
         , rank = 1, link = 'x---', dctype = 0                           &
         , dimension_ids = (/ DIMID_SEDI_PARTICLE /)                        &
         , ldimlen       = (/ nparticle_max /)                                &
         , axis = 'N---'                                                 &
         )
    CALL channel_halt(substr, status)


#ifdef MBM_CLAMS
   ! -----------------------------------------------------

    CALL new_representation(status, REPR_TIMEBNDS, 'TIME_BNDS'       &
         , rank = 1, link = 'x---', dctype = DC_BC                   &
         , dimension_ids = (/DIMID_BNDS/) &
         , ldimlen       = (/ 2    /)     &
         , axis = 'N---'                  &
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

    cnt(:,1)  = 2
    memu(:,1) = 2

    CALL set_representation_decomp(status, REPR_TIMEBNDS &
         , start, cnt, memu, meml, piotype=PIOTYPE_COL)
    CALL channel_halt(substr, status)    

    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

   ! -----------------------------------------------------
#endif    

    ! #################
    ! DIAGNOSTIC OUTPUT
    ! #################
    if (p_pe==0 .and. ldiagout) then
       CALL write_representation(status)
       CALL channel_halt(substr, status)
    endif

#ifdef MBM_CLAMS
  END SUBROUTINE main_channel_initialize_reprs
#else
  END SUBROUTINE clams_initialize_reprs
#endif
  ! -------------------------------------------------------------------
