!**********************************************************************
MODULE messy_clams_si

#if defined(ECHAM5) || defined(MBM_CLAMS)
!**********************************************************************
!  Submodel interface for CLaMS submodel 
!**********************************************************************
  ! SMCL
  USE messy_clams
  USE messy_clams_global,     ONLY: PREC, DP, &
                                    nmaxcltr  ! for clams gridding
  
  USE messy_main_timer_event, ONLY: time_event, io_time_event

!ju_ec_20180627+
! For CLaMS MESSy TRACER names
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM !TRACER BASENAME MAX LEN.
!ju_ec_20180627+

  IMPLICIT NONE
  SAVE ! op_pj_20180614

  ! MODULE VARIABLES
  TYPE(time_event)    :: clamsoutevent
  TYPE(io_time_event) :: io_clamsoutevent

  ! current airparcel positions and times
  REAL(PREC), DIMENSION(:), POINTER :: LAT        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LON        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LEV        => NULL()
  REAL(DP),   DIMENSION(:), POINTER :: JULSEC     => NULL()
  REAL(DP),   POINTER               :: JULTIME

  ! airparcel positions at last CHEM call
  REAL(PREC), DIMENSION(:), POINTER :: LAT_OLD        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LON_OLD        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LEV_OLD        => NULL()
  REAL(DP),   POINTER               :: JULTIME_OLD

  ! airparcel postions at last MIX/BMIX call
  REAL(PREC), DIMENSION(:), POINTER :: LAT_OLD_MIX    => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LON_OLD_MIX    => NULL()

  ! state of mixing / vertical mixing
  REAL(PREC), DIMENSION(:), POINTER :: STATE       => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: STATE_VERT  => NULL()

  ! for vertical mixing:
  REAL(PREC), DIMENSION(:), POINTER :: THETA_OLD_MIX  => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: BVF_OLD_MIX => NULL()
  ! THETA and BVF_WET must be specified on parameterlist (->PARAM),
  ! if vertical mixing is switched on

  ! Vertical grid used in MIX/BMIX:
  REAL(PREC), DIMENSION(:), POINTER :: LEV_GRID   => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LEV_DELTA  => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: R_GRID     => NULL()
 
  ! DRIVER VARIABLES (GLOBAL FIELD)
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_LAT3D_G  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_LON3D_G  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_UWIND_G  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_VWIND_G  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_WWIND_G  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_TEMP_G   => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_PRESS_G  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_THETA_G  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_OMEGA_G  => NULL()
  REAL(DP), DIMENSION(:,:)  , POINTER :: E5_PSURF_G  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_TTE_G    => NULL()
  REAL(DP), DIMENSION(:,:)  , POINTER :: E5_ALPSTE_G => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_THETADOT_G=> NULL() 
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_THETADOT2_G=> NULL() 
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_ZETA_G    => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_ZETADOT_G => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_SIGMA_G => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_SIGMADOT_G => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_HEATING_G => NULL()
  REAL(DP), DIMENSION(:,:)  , POINTER :: E5_PSDOT_G => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_SH_G => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_IWC_G => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_CLWC_G => NULL()
! op_pj_20180614+
!!$  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_HEAT_VDIFF_G => NULL()
!!$  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_HEAT_CLOUDCONV_G => NULL()
!!$  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_HEAT_RHEAT_G => NULL()
!!$  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_HEAT_SSO_G => NULL()
  ! diabatic temperature tendencies
#ifdef MESSYTENDENCY
  REAL(DP), DIMENSION(:,:,:), POINTER :: tte_dyn => NULL()
#endif
! op_pj_20180614-

!!!!! not used ???
!!$  REAL(DP), DIMENSION(:), POINTER :: hyam => NULL() ! Hybrid level A coefficients
!!$  REAL(DP), DIMENSION(:), POINTER :: hybm => NULL() ! Hybrid level B coefficients
!!$  REAL(DP), DIMENSION(:), POINTER :: hyai => NULL() ! Hybrid level A coefficients
!!$  REAL(DP), DIMENSION(:), POINTER :: hybi => NULL() ! Hybrid level B coefficients


  ! DRIVER VARIABLES (parallel decomposition)
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_LAT3D_D  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_LON3D_D  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_ZETA_D    => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_OMEGA_D  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_THETADOT_D=> NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_THETADOT2_D=> NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_ZETADOT_D => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_SIGMA_D => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_SIGMADOT_D => NULL()
  REAL(DP), DIMENSION(:,:)  , POINTER :: E5_PSDOT_D => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_HEATING_D => NULL()

  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_UWIND_D  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_VWIND_D  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_WWIND_D  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_TEMP_D   => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_PRESS_D  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_THETA_D  => NULL()
  REAL(DP), DIMENSION(:,:)  , POINTER :: E5_PSURF_D  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_TTE_D    => NULL()
  REAL(DP), DIMENSION(:,:)  , POINTER :: E5_ALPSTE_D => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_SH_D => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_IWC_D => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_CLWC_D => NULL()

! op_pj_20180614+
!!$  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_HEAT_VDIFF_D => NULL()
!!$  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_HEAT_CLOUDCONV_D => NULL()
!!$  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_HEAT_RHEAT_D => NULL()
!!$  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_HEAT_SSO_D => NULL()
! op_pj_20180614-
     
!!$  ! Pointer targets
!!$  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: UDT_T 
!!$  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: VDT_T 
!!$  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: WDT_T 

#ifdef ECHAM5
! ju_ec_20180627+
! CLaMS Tracer Set Information
  INTEGER, DIMENSION(NMAXCLTR)  :: CLTR_IDX !Indecies of CLaMS tracers
  INTEGER, DIMENSION(NMAXCLTR)  :: CLTR_SPECARR_IND !indexes of tracers 
                                                    !in specarr
  !
  ! The CLaMS tracers' indecies are defined in clams_init_coupling
  !
! ju_ec_20180627-

! ju_ec_20180717+
! Other additions for coupling to EMAC radiation
  INTEGER :: nlat , nlev, nlon  ! shape of ECHAM grid
  REAL(DP), DIMENSION(:), ALLOCATABLE :: E5_LAT_BOUNDS  ! Nr of lat bounds is 
                                                        ! one more than the nr 
                                                        ! of latitudes
  REAL(DP), DIMENSION(:), ALLOCATABLE :: E5_LON_BOUNDS  ! Nr of lon bounds is
                                                        ! same as nr of lons
! op_sb_20190712+
!!$  INTEGER, DIMENSION(:,:), ALLOCATABLE :: pos
!  REAL(DP), DIMENSION(:,:), POINTER :: POS   ! integer parcel positions moved to core
  INTEGER, PARAMETER:: POS_CART=3             ! 2. dimension for field POS(dnparts_max,:)
! op_sb_20190712-

!!  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: pc_d
!!  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: pc_g
! now channel objects
  REAL(DP), DIMENSION(:,:,:), POINTER :: pc_d => NULL()
!!  REAL(DP), DIMENSION(:,:,:), POINTER :: pc_g => NULL() ! moved to core
  ! pc_g in decomposition (channel object)
  REAL(DP), DIMENSION(:,:,:), POINTER :: spc_g => NULL()  ! op_sb_20190807
  ! op_sb_20190821
  REAL(DP), DIMENSION(:,:),   POINTER :: pblh_i  ! index of pblh from tropop

  REAL(DP), DIMENSION(:,:,:,:), POINTER :: gp_g => NULL()
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: gp_d => NULL()
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: clgp_data => NULL()

  ! op_sb_20200303+
  REAL(dp), DIMENSION(:,:,:), POINTER :: PRESSI_INIT  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: PRESSI_3D_G  => NULL()
  ! op_sb_20200303-

  ! ju_ec_20180717-
  ! op_sb_20191018+
    REAL(DP), DIMENSION(:,:,:), POINTER :: ptr_trac
    CHARACTER(LEN=40) :: name
    CHARACTER(LEN=40) :: longname
    CHARACTER(LEN=40) :: unit
  ! op_sb_20191018-
#endif


  PUBLIC :: clams_setup
  PUBLIC :: clams_initialize
  PUBLIC :: clams_init_memory
#ifdef ECHAM5
  PUBLIC :: clams_init_tracer ! op_sb_20200227
#endif
  PUBLIC :: clams_init_coupling
  PUBLIC :: clams_local_start    ! op_pj_20180614
  PUBLIC :: clams_local_end      ! op_pj_20180614
  PUBLIC :: clams_global_start
  PUBLIC :: clams_global_end
  PUBLIC :: clams_free_memory
#ifdef ECHAM5
! ju_ec_20180627+
  PUBLIC :: clams_new_tracer
! ju_ec_20180626-
#endif

  
!  PRIVATE

  PRIVATE :: calculate_thetadot_echam
  PRIVATE :: calculate_zeta_zetadot_echam

#ifdef ECHAM5
  PRIVATE :: clams_tracer_update                    ! ju_ec_20180629
  PRIVATE :: clams_lg2gp                            ! ju_ec_20180713
#endif  
!--------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clams_setup

    ! Call of nml reading subroutine, setup parallel decomposition

    ! Read some global information from init and first meteorological file
    
    ! Do some checks 

    USE messy_main_tools,          ONLY: find_next_free_unit
    USE messy_main_switch,         ONLY: USE_CLAMS, USE_CLAMSMIX, USE_CLAMSTRAJ, &
                                         USE_DISSOC, USE_CLAMSSEDI, &
                                         USE_CLAMSBMIX, USE_CLAMSCIRRUS, &
                                         USE_CLAMSDEEPCONV
!!#D clamschem +
    USE messy_main_switch,         ONLY: USE_CLAMSCHEM
!!#D clamschem -
    USE messy_main_timer,          ONLY: YEAR_START, MONTH_START, DAY_START, &
                                         MINUTE_START, HOUR_START, SECOND_START

    ! BMIL
    USE messy_main_blather_bi,     ONLY: error_bi
    USE messy_main_transform_bi,   ONLY: get_dc_index
    USE messy_main_mpi_bi,         ONLY: p_nprocs, p_pe, p_sum, p_barrier

    USE messy_main_timer,          ONLY: delta_time
    USE messy_clams_global,        ONLY: ntasks, rank, pi, lcoupled, resume_run, &
                                         nparts, nparts_max, dnparts, dnparts_max, rres, idx, &
                                         nparts_max_shuffle, dnparts_max_shuffle, rres_shuffle, &
                                         init_vertcoorname, met_vertcoorname, &
                                         initfile, first_initfile, ldiagout, &
                                         met_dir, met_prefix, theta_dir, theta_prefix, &
                                         pre_metfile, pre_thetafile, &
                                         pre_year, pre_month, pre_day, pre_sec, &
                                         asc_level, asc_lat, level_is_vertcoor, loglev, logpress, &
                                         nparams, nparams_old, paramnames, paramnames_old, &
                                         nx, ny, nz, ntheta, &
                                         latgrid, longrid, levelgrid, thetagrid, &
                                         latgrid_rad, longrid_rad, &
                                         corrfile, corr_thetadot, met_freq, buffersize
    USE messy_clams_tools_utils,   ONLY: lowercase, uppercase
    USE messy_clams_read_metdata,  ONLY: nc_read_corrfile
    USE messy_clams_tools_ncutils, ONLY: nc_get_vertcoorname, nc_get_level, nc_grid_descr
    USE messy_clams_tools_utils,   ONLY: uppercase, lowercase, str_found, str_pos,  &
                                         get_metfilename, bubble_sort

    USE netcdf

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr='clams_setup'
    INTEGER                     :: nlevels
    INTEGER                     :: status, iou, error, i, ipos, ipos2
    INTEGER                     :: rcode, ncid, dimid, ncid_first, varid
    INTEGER                     :: nparts_first
    INTEGER                     :: ntasks_init

    INTEGER, dimension(:), allocatable :: dnparts_init

    IF (p_pe==0) THEN
       WRITE(*,*)
       WRITE(*,*) uppercase(substr)
       WRITE(*,*)
    ENDIF

    ! Determine size of local arrays
    ntasks = p_nprocs
    rank = p_pe

#ifndef NOMPI
  ! original code with MPI
#else
  CALL error_bi('CLaMS cannot be run without MPI; please recompile with MPI',substr) 
#endif

    !--------------------------------------------------------------
    ! Coupled to ECHAM5 ?
    !--------------------------------------------------------------
#ifdef ECHAM5
    lcoupled = .true.
#else
    lcoupled = .false.
#endif

    !--------------------------------------------------------------
    ! Read namelist variables:
    !--------------------------------------------------------------

    iou = find_next_free_unit(100,200)
    CALL clams_read_nml(status, iou) 
    IF (status /= 0) CALL error_bi('Error in clams_read_nml !',substr)

    !--------------------------------------------------------------
    ! short names of output parameters
    !--------------------------------------------------------------
    IF (p_pe==0) THEN
       WRITE (*,*) 'Output parameters:'
       DO i = 1, nparams
          write (*,*) paramnames(i)
       ENDDO
       WRITE(*,*)
    ENDIF

    !--------------------------------------------------------------
    ! check parameters:
    !
    !  if CHEM is used:
    !      - TEMP and PRESS must be specified in paramlist
    !      - nparams_old is set to 2 (TEMP_OLD, PRESS_OLD)
    !
    !  if SEDI is used: TEMP and PRESS must be specified in paramlist
    !
    !  if EQLAT or PV is specified in paramlist: THETA must be specified too
    !
    !  if EQLAT or PV is specified in paramlist or BMIX is switched on: 
    !       theta_dir and theta_prefix must be specified too
    !
    !  if CIRRUS is used: 
    !      - TEMP and PRESS must be specified in clams.nml
    !      - PV must be specified in clams.nml
    !        (if PV is specified THETA must be specified too)
    !
    !  if MIX is used and vertical mixing is switched on:
    !       THETA and BVF_WET must be specified 
    !     => checked in clamsmix_initialize (after reading clamsmix.nml) !!!
    !
    ! if DEEPCONV is used:
    !    - TEMP, PRESS, THETA and BVF_WET  must be specified in clams.nml
    !
    !--------------------------------------------------------------

    !  if CHEM is used:
    !         TEMP and PRESS must be specified in paramlist
    !         nparams_old is set to 2 (TEMP_OLD, PRESS_OLD)
    nparams_old = 0 ! op_pj_20170110
!!#D clamschem +
    if (USE_CLAMSCHEM) then
       if (.NOT. str_found (nparams, paramnames, 'TEMP'))  &
            call error_bi ('Add parameter TEMP to paramlist (clams.nml) !!!',substr)
       if (.NOT. str_found (nparams, paramnames, 'PRESS')) &
            call error_bi ('Add parameter PRESS to paramlist (clams.nml) !!!',substr)
       nparams_old = 2
       paramnames_old(1) = 'TEMP'
       paramnames_old(2) = 'PRESS'
!!$    else                    ! op_pj_20170110
!!$       nparams_old = 0      ! op_pj_20170110
    endif
!!#D clamschem -

    !  if SEDI is used: TEMP and PRESS must be specified in paramlist
    if (USE_CLAMSSEDI) then
       if (.NOT. str_found (nparams, paramnames, 'TEMP'))  &
            call error_bi ('Add parameter TEMP to paramlist (clams.nml) !!!',substr)
       if (.NOT. str_found (nparams, paramnames, 'PRESS')) &
            call error_bi ('Add parameter PRESS to paramlist (clams.nml) !!!',substr)
    endif
   
#ifdef MBM_CLAMS    
    !  if EQLAT or PV is specified in paramlist: THETA must be specified too
    if (str_found (nparams, paramnames, 'EQLAT') .or. &
        str_found (nparams, paramnames, 'PV') ) then
       ipos = str_pos (nparams, paramnames, 'THETA')
       if (ipos<0) then
          call error_bi &
               ('If EQLAT or PV is used, THETA must be specified as parameter too (clams.nml) !!!', &
               substr)
       endif
    endif
    !  if EQLAT or PV is specified in paramlist:
    !       theta_dir and theta_prefix must be specified
    if (str_found (nparams, paramnames, 'EQLAT') .or. &
        str_found (nparams, paramnames, 'PV') ) then
       if (theta_dir=='' .or. theta_prefix=='') then
          call error_bi ('Specify theta_dir and theta_prefix used for EQLAT and PV (clams.nml) !!!',substr)
       endif
    endif
    !  if BMIX is switched on:
    !       theta_dir and theta_prefix must be specified
    if (USE_CLAMSBMIX) then
       if (theta_dir=='' .or. theta_prefix=='') then
          call error_bi &
               ('Specify theta_dir and theta_prefix used for BMIX (clams.nml) !!!',substr)
       endif
    endif   
    !  if CIRRUS is switch on: TEMP, PRESS, PV and THETA must be specified in clams.nml
    if (USE_CLAMSCIRRUS) then
       if (.NOT. str_found (nparams, paramnames, 'TEMP'))  &
            call error_bi ('Add parameter TEMP to paramlist (clams.nml) !!!',substr)
       if (.NOT. str_found (nparams, paramnames, 'PRESS')) &
            call error_bi ('Add parameter PRESS to paramlist (clams.nml) !!!',substr)
       ipos = str_pos (nparams, paramnames, 'PV')
       ipos2 = str_pos (nparams, paramnames, 'THETA')
       if (ipos<0 .or. ipos2<0) then
          call error_bi &
               ('If CIRRUS is switched on, PV and THETA must be specified as parameters in clams.nml !!!',substr)
       endif
    endif
#endif

    ! If DEEPCONV is switch on:
    ! TEMP, PRESS, THETA and BVF_WET must be specified in clams.nml
    if (USE_CLAMSDEEPCONV) then
       if (.NOT. str_found (nparams, paramnames, 'TEMP'))  &
            call error_bi ('DEEPCONV: Add parameter TEMP to paramlist (clams.nml) !!!',substr)
       if (.NOT. str_found (nparams, paramnames, 'PRESS'))  &
            call error_bi ('DEEPCONV: Add parameter PRESS to paramlist (clams.nml) !!!',substr)
       if (.NOT. str_found (nparams, paramnames, 'THETA'))  &
            call error_bi ('DEEPCONV: Add parameter THETA to paramlist (clams.nml) !!!',substr)
       if (.NOT. str_found (nparams, paramnames, 'BVF_WET'))  &
            call error_bi ('DEEPCONV: Add parameter BVF_WET to paramlist (clams.nml) !!!',substr)
    endif
    
    !--------------------------------------------------------------
    ! check module dependencies:
    !
    ! - CLAMS must be switched on
    ! - CLAMSTRAJ must be switched on (?)
    ! - If CLAMSCHEM is switched on, DISSOC must be switched on too
    ! - If CLAMSSEDI is switched on, CLAMSCHEM must be switched on too
    ! - If CLAMSBMIX is switched on, CLAMSMIX must be switched on too
    ! - if CLAMSDEEPCONV is switched on, CLAMSCHEM must be switched on too
    !  
    ! Note:
    ! If CLAMSMIX is switched on, mixing can be disabled with 
    ! switch_mixing=0 !!!
    !
    !--------------------------------------------------------------
    if (.not. USE_CLAMS) then
       call error_bi ('CLAMS must be switched on !!!',substr)
    endif
    if (.not. USE_CLAMSTRAJ) then
       call error_bi ('CLAMSTRAJ must be switched on !!!',substr)
    endif
!!#D clamschem +
    if (USE_CLAMSCHEM) then
       if (.not. USE_DISSOC) &
            call error_bi ('If CHEM is switched on, DISSOC must be switched on too !!!',substr)
    endif
    if (USE_CLAMSSEDI) then
       if (.not. USE_CLAMSCHEM) &
            call error_bi ('If SEDI is switched on, CHEM must be switched on too !!!',substr)
    endif
!!#D clamschem -
    if (USE_CLAMSBMIX) then
       if (.not. USE_CLAMSMIX) &
            call error_bi ('If BMIX is switched on, MIX must be switched on too !!!',substr)
    endif
    if (USE_CLAMSDEEPCONV)then
       if (.not. USE_CLAMSCHEM) &
            call error_bi ('If DEEPCONV is switched on, CHEM must be switched on too !!!',substr)
    endif

   
    !--------------------------------------------------------------
    ! check timesteps
    !--------------------------------------------------------------
    if (mod(met_freq*3600,int(delta_time)) /= 0) then
       call error_bi ("Interval between windfiles (met_freq) must be multiple of delta_time !!!",substr)
    endif
   


    !--------------------------------------------------------------
    ! Read from init file:
    !
    ! - number of airparcels (NPARTS)
    ! - Name of vertical coordinate
    ! - number of theta/zeta-levels
    !
    ! For resume of run (NOT messy restart!) with new initfile:
    ! - ntasks
    ! - dnparts(ntasks)
    !--------------------------------------------------------------

    if (first_initfile == '') first_initfile = initfile

    
    ! Open init files
    rcode = nf90_open (initfile, nf90_nowrite, ncid, buffersize)
    IF (rcode /= 0) CALL error_bi('Cannot open file '//trim(initfile),substr)
    rcode = nf90_open (first_initfile, nf90_nowrite, ncid_first, buffersize)
    IF (rcode /= 0) CALL error_bi('Cannot open file '//trim(first_initfile),substr)

    ! Read nparts from init file
    if (p_pe==0) write (*,*) 'Read nparts from ',trim(initfile)
    rcode = nf90_inq_dimid (ncid,"NPARTS",dimid)
    IF (rcode /= 0) CALL error_bi('Cannot find dimension nparts !',substr)
    rcode = nf90_inquire_dimension (ncid,dimid,len=nparts)
    IF (rcode /= 0) CALL error_bi('Cannot read dimension nparts !',substr)
    if (p_pe==0) write (*,*) 'Read nparts_first from ',trim(first_initfile)
    rcode = nf90_inq_dimid (ncid_first,"NPARTS",dimid)
    IF (rcode /= 0) CALL error_bi('Cannot find dimension nparts !',substr)
    rcode = nf90_inquire_dimension (ncid_first,dimid,len=nparts_first)
    IF (rcode /= 0) CALL error_bi('Cannot read dimension nparts !',substr)
    if (p_pe==0) WRITE (*,*) 'NPARTS=',nparts
    if (p_pe==0) WRITE (*,*) 'NPARTS_FIRST=',nparts_first
    if (p_pe==0) WRITE (*,*)
    

    ! Name of vertical coordinate in init file
    if (p_pe==0) write (*,*) 'Read name of vertical coordinate from ',trim(initfile)
    if (p_pe==0) WRITE (*,*)
    call nc_get_vertcoorname (ncid, init_vertcoorname)

    ! Read number of levels (NTHETAS/NZETAS) 
    if ((trim(lowercase(init_vertcoorname))=="theta" .or. &
         trim(lowercase(init_vertcoorname))=="zeta") .and. USE_CLAMSMIX) then
       rcode = nf90_inq_dimid(ncid,"N"//trim(uppercase(init_vertcoorname))//"S",dimid)
       IF (rcode /= 0) CALL error_bi &
            ('Cannot find dimension '//"N"//trim(uppercase(init_vertcoorname))//"S",substr)
       rcode = nf90_inquire_dimension (ncid,dimid,len=nlevels)
       IF (rcode /= 0) CALL error_bi &
            ('Cannot read dimension '//"N"//trim(uppercase(init_vertcoorname))//"S",substr)
    else if (trim(lowercase(init_vertcoorname))=="press" .and. USE_CLAMSMIX)&
         & then
       !!! Option fuer NTHETAS ergaenzen!!!!
       rcode = nf90_inq_dimid(ncid,"NZETAS",dimid)
       IF (rcode /= 0) CALL error_bi ('Cannot find dimension NZETAS',substr)
       rcode = nf90_inquire_dimension (ncid,dimid,len=nlevels)
       IF (rcode /= 0) CALL error_bi ('Cannot read dimension NZETAS',substr)
    else
       nlevels = -1
    endif

    ! resume run with new init file: read dnparts(ntasks)
    if (resume_run) then

       ! read and check dimension ntasks 
       rcode = nf90_inq_dimid (ncid,"ntasks",dimid)
       if (rcode /= 0) call error_bi('Cannot find dimension ntasks in init file!',substr)
       rcode = nf90_inquire_dimension (ncid,dimid,len=ntasks_init)
       if (rcode /= 0) call error_bi('Cannot read dimension ntasks from init file!',substr)
       if (ntasks_init /= ntasks) call error_bi('NTASKS has changed  !',substr)
       
       ! read dnparts(ntasks)
       rcode = nf90_inq_varid (ncid,"dnparts",varid)
       if (rcode /= 0) call error_bi('Cannot find variable dnparts in init file!',substr)
       allocate (dnparts_init(ntasks))
       rcode = nf90_get_var (ncid, varid, dnparts_init)
       if (rcode /= 0) call error_bi('Cannot read variable dnparts from init file!',substr)

       if(p_pe==0) write (*,*) 'dnparts_init=',dnparts_init

    endif
       
    ! Close init files
    rcode = nf90_close (ncid)
    IF (rcode /= 0) CALL error_bi('Cannot close file '//trim(initfile),substr)
    rcode = nf90_close (ncid_first)
    IF (rcode /= 0) CALL error_bi('Cannot close file '//trim(first_initfile),substr)


    !--------------------------------------------------------------
    ! Set dimensions:
    !
    ! - nparts_max  : max. number of airparcels
    ! - nparts      : current number of airparcels (read from init file)
    ! - dnparts_max : max. number of airparcels on rank
    ! - dnparts     : current number of airparcels on rank
    !--------------------------------------------------------------

    ! Get number of cells per pe
    CALL get_dc_index(nparts, idx)
       
    dnparts  = idx(p_pe,2) - idx(p_pe,1) + 1
       
    if (p_pe==0) write (*,*) 'resume_run=', resume_run

    if (resume_run) then

       ! p_sum(dnparts) = nparts ?

       idx(0,1) = 1
       idx(0,2) = dnparts_init(1)
       do i = 1, p_nprocs-1
          idx(i,1) = idx(i-1,2)+1
          idx(i,2) = idx(i-1,2)+dnparts_init(i+1)
       enddo

       dnparts = dnparts_init(p_pe+1)
       write (*,*) 'rank, dnparts=',p_pe,dnparts

       deallocate (dnparts_init)

    endif

    if (nlevels == -1) then
       dnparts_max = nparts_first / ntasks * rres
       dnparts_max_shuffle = nparts_first / ntasks * rres_shuffle
    else
       dnparts_max = nparts_first / min(ntasks,nlevels) * rres
       dnparts_max_shuffle = nparts_first / min(ntasks,nlevels) * rres_shuffle
    endif
    
    nparts_max  = p_sum(dnparts_max)
    nparts_max_shuffle  = p_sum(dnparts_max_shuffle)
    
    if (dnparts > dnparts_max) then
       if (p_pe==0) then
          write(*,*) 'dnparts     ', dnparts
          write(*,*) 'dnparts_max ', dnparts_max
       endif
       CALL error_bi('dnparts > dnparts_max !',substr)
    endif
    
    if (p_pe==0) then
       write(*,*) 'p_pe', p_pe, 'dnparts     ', dnparts
       write(*,*) 'p_pe', p_pe, 'dnparts_max ', dnparts_max
       write(*,*) 'p_pe', p_pe, 'nparts      ', nparts
       write(*,*) 'p_pe', p_pe, 'nparts_max  ', nparts_max
       write(*,*) 'p_pe', p_pe, 'nparts_max_shuffle  ', nparts_max_shuffle
       write(*,*) 'p_pe', p_pe, 'dnparts_max_shuffle ', dnparts_max_shuffle
    endif

    
#ifdef MBM_CLAMS    

    !--------------------------------------------------------------
    ! Read grid from first meteorological file:
    ! - nx, ny, nz
    ! - latgrid, longrid, levelgrid
    ! Set logical variables:
    ! - level_is_vertcoor
    ! - asc_level
    ! - loglev
    ! - logpress
    !--------------------------------------------------------------

!!!!! Funktioniert nur, wenn Startzeit=Datenzeit

    pre_year  = YEAR_START
    pre_month = MONTH_START
    pre_day   = DAY_START
    pre_sec   = HOUR_START*3600 + MINUTE_START*60 + SECOND_START 

    ! get name of first meteorological file
    pre_metfile = get_metfilename (met_prefix, met_dir, &
                       pre_year, pre_month, pre_day, pre_sec/3600)
    if (p_pe==0) write(*,*) 'pre_metfile', pre_metfile
    if (p_pe==0) write(*,*) 'read vertical coordinate from met_file', pre_metfile
    if (theta_prefix=='') then
       pre_thetafile = ''
    else
       pre_thetafile = get_metfilename (theta_prefix, theta_dir, &
            pre_year, pre_month, pre_day, pre_sec/3600)
       if (p_pe==0) write(*,*) 'pre_thetafile', pre_thetafile
    endif

    ! Name of vertical coordinate in meteorological datasets
    call nc_get_vertcoorname (pre_metfile, met_vertcoorname)
    met_vertcoorname = ADJUSTL(met_vertcoorname)


    ! check, if vertical coordinate in init-file and windfile are different
    if (TRIM(lowercase(met_vertcoorname)) == TRIM(lowercase(init_vertcoorname))) then
       level_is_vertcoor = .true.
    else
       level_is_vertcoor = .false.
    endif
  
    ! interpolation for levels linear or log. linear 
    loglev = .false.     ! default: linear
    asc_level = .true.   ! default: ascending levels
    
    ! if level (vertical coordinate in pos-files) is PRESS:
    if (trim(lowercase(init_vertcoorname)) == 'press') then
       loglev = .true.      ! log. linear
       asc_level = .false.  ! decending levels
    endif
     
    IF(p_pe==0)THEN
       write (*,*)
       write (*,*) 'Use vertical coordinate ',trim(init_vertcoorname)
       write (*,*) 
       write (*,*) 'Vertical coordinate in windfile: ',trim(met_vertcoorname)
       write (*,*)
       write (*,*) 'Interpolation for levels logarithmic: ', loglev
       write (*,*)
    ENDIF

!!!!! Hier aendern, wenn PRESS linear interpoliert werden soll:
    logpress = .true.

    if(p_pe==0)THEN
       if (trim(init_vertcoorname) /= 'press') then
          if (logpress) then
             write (*,*)
             write (*,*) 'Interpolation for PRESS:  '
             write (*,*) '     <=500K:  linear'
             write (*,*) '     >=1000K: logarithmic'
             write (*,*) '     >500K and <1000K: linear/logarithmic'
             write (*,*)
          else
             write (*,*)
             write (*,*) 'Linear Interpolation for PRESS'
             write (*,*)
          endif
       endif
    endif
    
    ! Determine levels 
    if (p_pe==0) WRITE (*,*) 'Read level from file ', TRIM(pre_metfile)
    CALL nc_get_level (TRIM(pre_metfile), nz, levelgrid, err=error)
    IF (error /= 0 ) call error_bi ("Level could not be read !!!",substr)
    if (p_pe==0) WRITE (*,*) 'nz=',nz
    if (pre_thetafile == '') then
       ntheta = 0
    else
       if (p_pe==0) WRITE (*,*) 'Read level from file ', TRIM(pre_thetafile)
       CALL nc_get_level (trim(pre_thetafile), ntheta, thetagrid, err=error)
       IF (error /= 0 ) call error_bi ("Theta-level could not be read !!!",substr)
       if (p_pe==0) WRITE (*,*) 'ntheta=',ntheta
    endif
    
 
    ! Extract grid description of wind data
    CALL nc_grid_descr(pre_metfile,nx,ny,longrid,latgrid,status) 

    CALL bubble_sort (levelgrid,nz)

    if (p_pe==0) then
       WRITE (*,*)
       WRITE (*,'(A,I4)') 'Number of longitudes: ', nx
       WRITE (*,'(A,I4)') 'Number of latitudes: ', ny
       WRITE (*,'(A,I4)') 'Number of levels: ', nz
       WRITE (*,*)
       if (ldiagout) then
          write (*,*) 'longrid=',longrid
          write (*,*) 'latgrid=',latgrid
          write (*,*) 'levelgrid=',levelgrid
          WRITE (*,*)
       endif
    endif
 
    if (longrid(1) < 0.) call error_bi('longitudes are not in valid range !!!',substr)

    allocate (longrid_rad(nx+1))
    allocate (latgrid_rad(0:ny+1))

    longrid_rad(1:nx) = longrid(1:nx)/180. * pi
    latgrid_rad(1:ny) = latgrid(1:ny)/180. * pi
    
    longrid_rad(nx+1) = (longrid(1)+360.)/180. * pi
    if (latgrid(2)>latgrid(1)) then
       asc_lat = .true.
       latgrid_rad(0)    = -0.5 * pi
       latgrid_rad(ny+1) = 0.5 * pi
    else
       asc_lat = .false.
       latgrid_rad(0)    = 0.5 * pi
       latgrid_rad(ny+1) = -0.5 * pi
    endif

#endif

#ifdef ECHAM5
    
    !--------------------------------------------------------------
    ! Set logical variables:
    ! - level_is_vertcoor
    ! - asc_level
    ! - loglev
    ! - logpress
    !--------------------------------------------------------------
    level_is_vertcoor = .false.   ! model levels are used !
    asc_level = .true. 
    loglev = .false.
    logpress = .true.
    
!!!!! nx, ny, nlev, latgrid, longrid, levelgrid, longrid_rad, latgrid_rad 
!!!!! => in clams_global_end ?!?  

#endif    


!!!!! ???
#ifdef MBM_CLAMS    
    ! If DEEPCONV is switched on:
    ! ZETA must be vertical coordinate in initfile
    ! windfiles on hybrid levels
    if (USE_CLAMSDEEPCONV) then
       if (trim(uppercase(init_vertcoorname)) /= 'ZETA') &
            call error_bi('DEEPCONV: ZETA must be vertical coordinate in initfile !!!',substr)
       if (trim(uppercase(met_vertcoorname)) /= 'HYBRID') &
            call error_bi('DEEPCONV: vertical coordinate in windfiles must be HYBRID!',substr)
    endif
#endif    
           


#ifdef MBM_CLAMS
   !--------------------------------------------------------------
   ! correction file for thetadot
   !--------------------------------------------------------------
   if (corrfile=='') then
      if (p_pe==0) write (*,*) 'no correction file'
      if (p_pe==0) write (*,*) 
      corr_thetadot = .false.
   else
      if (p_pe==0) write (*,*) 'correction file: ',trim(corrfile)
      if (p_pe==0) write (*,*) 
      call nc_read_corrfile (status,corrfile,init_vertcoorname)
      IF (status /= 0) CALL error_bi('Cannot read correction file !',substr)
      corr_thetadot = .true.
   endif
#endif

  END SUBROUTINE clams_setup

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clams_initialize

    ! BMIL
    USE messy_main_mpi_bi,      ONLY: p_pe
    USE messy_main_blather_bi,  ONLY: error_bi
    USE messy_main_timer_bi,    ONLY: timer_event_init

    ! SMCL
    USE messy_main_timer,        ONLY: delta_time
    USE messy_clams_global,      ONLY: SPECARR, maxspec, specnames, &
                                       PARAM, PARAM_OLD, nparams, nparams_old, &
                                       paramnames, paramnames_old, &
                                       PREDATA, FUTDATA, &
                                       timestep_clamsout,  &
                                       dnparts_max   ! op_sb_209190716
    USE messy_clams_tools_utils, ONLY: uppercase
#ifdef ECHAM5    
 ! op_sb_20190716+
    USE messy_main_channel_bi,         ONLY: DIMID_POS, DIMID_TRAJ, REPR_POS, DC_IX
    USE messy_main_channel_error_bi,   ONLY: channel_halt       
    USE messy_main_channel_dimensions, ONLY: new_dimension
    USE messy_main_channel_repr,       ONLY: new_representation, AUTO 
 ! op_sb_20190716-
#endif

    IMPLICIT NONE

    
    INTEGER :: status   ! op_sb_20190716
    CHARACTER(LEN=*), PARAMETER :: substr = 'clams_initialize'
    INTEGER :: i

    IF (p_pe==0) THEN
       WRITE(*,*)
       WRITE(*,*) uppercase(substr)
       WRITE(*,*)
    ENDIF

    ALLOCATE (SPECARR(MAXSPEC))
    ALLOCATE (specnames(MAXSPEC)) 

    IF (nparams > 0) THEN

       allocate (PARAM(nparams))
       allocate (PREDATA(nparams))
       allocate (FUTDATA(nparams))

       DO i = 1, nparams
          if (p_pe==0) write (*,*) 'i, paramnames(i)=',i,paramnames(i)
          PARAM(i)%name    = paramnames(i)
          PREDATA(i)%name  = paramnames(i)
          FUTDATA(i)%name  = paramnames(i)
       ENDDO
    ENDIF

    IF (nparams_old > 0) then
       allocate (PARAM_OLD(nparams_old))
       DO i = 1, nparams_old
          PARAM_OLD(i)%name = trim(paramnames_old(i))//'_OLD'
       ENDDO
    ENDIF
  
!!!!!
    if (mod(timestep_clamsout*3600,int(delta_time)) /= 0) then
       call error_bi ("CLAMS output timestep must be multiple of delta_time !!!",substr)
    endif

    ! Define CLAMS event:
    io_clamsoutevent%counter = timestep_clamsout
    io_clamsoutevent%unit = 'hours'
    io_clamsoutevent%adjustment = 'exact'
    io_clamsoutevent%offset = -delta_time
    CALL timer_event_init (clamsoutevent, io_clamsoutevent, 'CLAMSOUT_Event', 'present')

#ifdef ECHAM5
    CALL clams_initialize_gatts
    CALL clams_initialize_dims
    CALL clams_initialize_reprs
    CALL new_dimension(status, DIMID_POS, 'pos_cart', pos_cart)
    CALL channel_halt(substr, status)
    CALL new_representation(status, REPR_POS, 'REPR_POS'   &
         , rank = 2, link = 'xx--', dctype = DC_IX         &
         , dimension_ids = (/ DIMID_TRAJ, DIMID_POS /)     &
         , ldimlen       = (/ dnparts_max, AUTO /)         &
         )
    CALL channel_halt(substr, status)
#endif
  END SUBROUTINE clams_initialize
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clams_init_memory

    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_channel_error_bi, ONLY: channel_halt
#ifdef MBM_CLAMS    
    USE messy_main_channel_bi, ONLY: REPR_LG_CLAMS, &
                                     REPR_3DINP_CLAMS,  &
                                     REPR_3DINP_CLAMSTHETA, REPR_MIX_GRID, &
                                     REPR_NTASKS, SCALAR
#endif
#ifdef ECHAM5
    USE messy_main_channel_bi, ONLY: REPR_LG_CLAMS, &
                                     REPR_3DINP_CLAMS,  &
                                     REPR_3DINP_CLAMSTHETA, REPR_MIX_GRID, &
                                     REPR_NTASKS, SCALAR, &
                                     GP_3D_MID, GP_2D_HORIZONTAL,&
                                     REPR_POS   ! op_sb_20190716
#endif    
    USE messy_main_blather_bi, ONLY: error_bi

    ! SMCL
    USE messy_main_channel,    ONLY: new_channel, new_attribute, &
                                     new_channel_object, &
                                     get_channel_object_dimvar, & ! ju_ec_20180717
                                     new_channel_object_reference, & ! op_sb_20191021
                                     get_channel_info ! op_sb_20191022
    USE messy_main_switch,     ONLY: USE_CLAMSSEDI, USE_CLAMSMIX, USE_CLAMSBMIX, &
                                     USE_CLAMSDEEPCONV
    USE messy_clams,           ONLY: modstr, modver
    USE messy_clams_global,    ONLY: SPECARR, nspec, &
                                     PARAM, PARAM_OLD, nparams, nparams_old, &
                                     UDT, VDT, WDT, LEVELDT, DLEVDZDT, &
                                     UFUT, VFUT, WFUT, LEVELFUT, DLEVDZFUT, &
                                     PREDATA, FUTDATA, &
                                     username, mdi, &
                                     pre_metfile, pre_thetafile, &
                                     dnparts_co, dnparts_max_co, grid_switch_co, &
                                     pre_year_co, pre_month_co, &
                                     pre_day_co, pre_sec_co, &
                                     ldiagout, DP, sample_interval, init_vertcoorname,&
                                     E5_LAT, E5_LON,&       ! ju_ec_20180717  
                                     E5_LEVEL, dnparts_max,&! ju_ec_20180717
                                     n_cltr, &
                                     clams_gridding, clams_grid_verbose, rank, &
                                     MAXSPEC
    USE messy_clamsmix_global,     ONLY: switch_mixing, adapt_par
    USE messy_clamssedi_global,    ONLY: UDT_sedi, VDT_sedi, WDT_sedi, &
                                         leveldt_sedi, dlevdzdt_sedi
    USE messy_clams_tools_utils,   ONLY: uppercase
    USE messy_clams_tools_ncutils, ONLY: nc_get_var_atts
    USE messy_main_tools,          ONLY: PTR_1D_ARRAY ! ju_ec_20180717
    USE messy_main_constants_mem,  ONLY: STRLEN_ULONG ! ju_ec_20180717
    USE messy_main_tracer_mem_bi,  ONLY: ntrac_cl ! op_sb_20191018
    USE messy_main_tracer,         ONLY: get_tracer ! op_sb_20191018

    IMPLICIT NONE

    ! LOCAL 
    CHARACTER(LEN=*), PARAMETER :: substr = 'clams_init_memory'
    INTEGER       :: status, ispec, i, posloc, jt
    CHARACTER(40) :: cr_date
    CHARACTER(8)  :: ydate
    CHARACTER(10) :: ytime
    REAL(DP)      :: valid_min, valid_max


    TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER :: dva  => NULL()
    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER :: units => NULL()

    IF (p_pe==0) WRITE(*,*) uppercase(substr)

    CALL DATE_AND_TIME(ydate, ytime)
    WRITE(cr_date,'(A4,5(A,A2))')   &
         ydate(1:4),'-',ydate(5:6),'-',ydate(7:8),' ',  &
         ytime(1:2),':',ytime(3:4),':',ytime(5:6) 

    !-----------------------------------------------------------------
    ! Define channel CLAMS  
    !-----------------------------------------------------------------

    CALL new_channel  (status, modstr, lrestreq=.TRUE.)
    CALL channel_halt (substr, status)
    CALL new_attribute(status, modstr, modstr//'_version', c = modver)
    CALL channel_halt (substr, status)

    ! Define channel object for dnparts
    CALL new_channel_object(status, modstr, "dnparts", p1=dnparts_co, &
                            reprid=REPR_NTASKS, lrestreq=.TRUE.)
    CALL channel_halt (substr, status)
    ! Define channel object for dnparts_max
    CALL new_channel_object(status, modstr, "dnparts_max", p1=dnparts_max_co, &
                            reprid=REPR_NTASKS, lrestreq=.TRUE.)
    CALL channel_halt (substr, status)

    ! Define channel object for grid_switch
    CALL new_channel_object(status, modstr, "grid_switch", p0=grid_switch_co, &
                            reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt (substr, status)

    ! Define channel objects for next uvfile date
    CALL new_channel_object(status, modstr, "pre_year", p0=pre_year_co, &
                            reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt (substr, status)
    CALL new_channel_object(status, modstr, "pre_month", p0=pre_month_co, &
                            reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt (substr, status)
    CALL new_channel_object(status, modstr, "pre_day", p0=pre_day_co, &
                            reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt (substr, status)
    CALL new_channel_object(status, modstr, "pre_sec", p0=pre_sec_co, &
                            reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt (substr, status)



    CALL new_channel_object(status, modstr, 'JULTIME', p0=JULTIME, reprid=SCALAR)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'JULTIME' &
         , 'creator_of_parameter', c = TRIM(username))
    CALL new_attribute(status, modstr, 'JULTIME' &
         , 'param_creation_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'JULTIME' &
         , 'param_modification_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'JULTIME' &
         , 'long_name', c = 'Time')
    CALL new_attribute(status, modstr, 'JULTIME' &
         , 'units', c = 'seconds since 2000-01-01 00:00:00 UTC')
    CALL new_attribute(status, modstr, 'JULTIME' &
         , 'missing_value', r = mdi)
    CALL new_attribute(status, modstr, 'JULTIME' &
         , 'sample_interval', r = sample_interval)
    CALL new_attribute(status, modstr, 'JULTIME' &
         , 'flag', c = 'NONE')
    CALL new_attribute(status, modstr, 'JULTIME' &
         , 'description', c = 'Time in julian seconds')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'JULTIME_OLD', p0=JULTIME_OLD, reprid=SCALAR)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'JULTIME_OLD' &
         , 'creator_of_parameter', c = TRIM(username))
    CALL new_attribute(status, modstr, 'JULTIME_OLD' &
         , 'param_creation_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'JULTIME_OLD' &
         , 'param_modification_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'JULTIME_OLD' &
         , 'long_name', c = 'Time')
    CALL new_attribute(status, modstr, 'JULTIME_OLD' &
         , 'units', c = 'seconds since 2000-01-01 00:00:00 UTC')
    CALL new_attribute(status, modstr, 'JULTIME_OLD' &
         , 'missing_value', r = mdi)
    CALL new_attribute(status, modstr, 'JULTIME_OLD' &
         , 'sample_interval', r = sample_interval)
    CALL new_attribute(status, modstr, 'JULTIME_OLD' &
         , 'flag', c = 'NONE')
    CALL new_attribute(status, modstr, 'JULTIME_OLD' &
         , 'description', c = 'Time in julian seconds')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'JULSEC', p1=JULSEC, reprid=REPR_LG_CLAMS)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'JULSEC' &
         , 'creator_of_parameter', c = TRIM(username))
    CALL new_attribute(status, modstr, 'JULSEC' &
         , 'param_creation_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'JULSEC' &
         , 'param_modification_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'JULSEC' &
         , 'long_name', c = 'Time')
    CALL new_attribute(status, modstr, 'JULSEC' &
         , 'units', c = 'seconds since 2000-01-01 00:00:00 UTC')
    CALL new_attribute(status, modstr, 'JULSEC' &
         , 'missing_value', r = mdi)
    CALL new_attribute(status, modstr, 'JULSEC' &
         , 'sample_interval', r = sample_interval)
    CALL new_attribute(status, modstr, 'JULSEC' &
         , 'flag', c = 'NONE')
    CALL new_attribute(status, modstr, 'JULSEC' &
         , 'description', c = 'Time in julian seconds')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'LAT', p1=LAT, reprid=REPR_LG_CLAMS)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'LAT' &
         , 'creator_of_parameter', c = TRIM(username))
    CALL new_attribute(status, modstr, 'LAT' &
         , 'param_creation_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'LAT' &
         , 'param_modification_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'LAT' &
         , 'long_name', c = 'Latitude')
    CALL new_attribute(status, modstr, 'LAT' &
         , 'units', c = 'deg N')
    CALL new_attribute(status, modstr, 'LAT' &
         , 'missing_value', r = mdi)
    CALL new_attribute(status, modstr, 'LAT' &
         , 'sample_interval', r = sample_interval)
    CALL new_attribute(status, modstr, 'LAT' &
         , 'flag', c = 'NONE')
    valid_min = -90.
    valid_max = 90.
    CALL new_attribute(status, modstr, 'LAT' &
         , 'valid_min', r = valid_min)
    CALL new_attribute(status, modstr, 'LAT' &
         , 'valid_max', r = valid_max)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'LAT_OLD', p1=LAT_OLD, reprid=REPR_LG_CLAMS)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'LAT_OLD' &
         , 'creator_of_parameter', c = TRIM(username))
    CALL new_attribute(status, modstr, 'LAT_OLD' &
         , 'param_creation_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'LAT_OLD' &
         , 'param_modification_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'LAT_OLD' &
         , 'long_name', c = 'Latitude')
    CALL new_attribute(status, modstr, 'LAT_OLD' &
         , 'units', c = 'deg N')
    CALL new_attribute(status, modstr, 'LAT_OLD' &
         , 'missing_value', r = mdi)
    CALL new_attribute(status, modstr, 'LAT_OLD' &
         , 'sample_interval', r = sample_interval)
    CALL new_attribute(status, modstr, 'LAT_OLD' &
         , 'flag', c = 'NONE')
    valid_min = -90.
    valid_max = 90.
    CALL new_attribute(status, modstr, 'LAT_OLD' &
         , 'valid_min', r = valid_min)
    CALL new_attribute(status, modstr, 'LAT_OLD' &
         , 'valid_max', r = valid_max)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'LAT_OLD_MIX', p1=LAT_OLD_MIX, reprid=REPR_LG_CLAMS)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'LAT_OLD_MIX' &
         , 'creator_of_parameter', c = TRIM(username))
    CALL new_attribute(status, modstr, 'LAT_OLD_MIX' &
         , 'param_creation_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'LAT_OLD_MIX' &
         , 'param_modification_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'LAT_OLD_MIX' &
         , 'long_name', c = 'Latitude')
    CALL new_attribute(status, modstr, 'LAT_OLD_MIX' &
         , 'units', c = 'deg N')
    CALL new_attribute(status, modstr, 'LAT_OLD_MIX' &
         , 'missing_value', r = mdi)
    CALL new_attribute(status, modstr, 'LAT_OLD_MIX' &
         , 'sample_interval', r = sample_interval)
    CALL new_attribute(status, modstr, 'LAT_OLD_MIX' &
         , 'flag', c = 'NONE')
    valid_min = -90.
    valid_max = 90.
    CALL new_attribute(status, modstr, 'LAT_OLD_MIX' &
         , 'valid_min', r = valid_min)
    CALL new_attribute(status, modstr, 'LAT_OLD_MIX' &
         , 'valid_max', r = valid_max)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'LON', p1=LON, reprid=REPR_LG_CLAMS)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'LON' &
         , 'creator_of_parameter', c = TRIM(username))
    CALL new_attribute(status, modstr, 'LON' &
         , 'param_creation_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'LON' &
         , 'param_modification_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'LON' &
         , 'long_name', c = 'Longitude')
    CALL new_attribute(status, modstr, 'LON' &
         , 'units', c = 'deg E')
    CALL new_attribute(status, modstr, 'LON' &
         , 'missing_value', r = mdi)
    CALL new_attribute(status, modstr, 'LON' &
         , 'sample_interval', r = sample_interval)
    CALL new_attribute(status, modstr, 'LON' &
         , 'flag', c = 'NONE')
    valid_min = 0.
    valid_max = 360.
    CALL new_attribute(status, modstr, 'LON' &
         , 'valid_min', r = valid_min)
    CALL new_attribute(status, modstr, 'LON' &
         , 'valid_max', r = valid_max)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'LON_OLD', p1=LON_OLD, reprid=REPR_LG_CLAMS)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'LON_OLD' &
         , 'creator_of_parameter', c = TRIM(username))
    CALL new_attribute(status, modstr, 'LON_OLD' &
         , 'param_creation_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'LON_OLD' &
         , 'param_modification_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'LON_OLD' &
         , 'long_name', c = 'Longitude')
    CALL new_attribute(status, modstr, 'LON_OLD' &
         , 'units', c = 'deg E')
    CALL new_attribute(status, modstr, 'LON_OLD' &
         , 'missing_value', r = mdi)
    CALL new_attribute(status, modstr, 'LON_OLD' &
         , 'sample_interval', r = sample_interval)
    CALL new_attribute(status, modstr, 'LON_OLD' &
         , 'flag', c = 'NONE')
    valid_min = 0.
    valid_max = 360.
    CALL new_attribute(status, modstr, 'LON_OLD' &
         , 'valid_min', r = valid_min)
    CALL new_attribute(status, modstr, 'LON_OLD' &
         , 'valid_max', r = valid_max)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'LON_OLD_MIX', p1=LON_OLD_MIX, reprid=REPR_LG_CLAMS)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'LON_OLD_MIX' &
         , 'creator_of_parameter', c = TRIM(username))
    CALL new_attribute(status, modstr, 'LON_OLD_MIX' &
         , 'param_creation_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'LON_OLD_MIX' &
         , 'param_modification_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'LON_OLD_MIX' &
         , 'long_name', c = 'Longitude')
    CALL new_attribute(status, modstr, 'LON_OLD_MIX' &
         , 'units', c = 'deg E')
    CALL new_attribute(status, modstr, 'LON_OLD_MIX' &
         , 'missing_value', r = mdi)
    CALL new_attribute(status, modstr, 'LON_OLD_MIX' &
         , 'sample_interval', r = sample_interval)
    CALL new_attribute(status, modstr, 'LON_OLD_MIX' &
         , 'flag', c = 'NONE')
    valid_min = 0.
    valid_max = 360.
    CALL new_attribute(status, modstr, 'LON_OLD_MIX' &
         , 'valid_min', r = valid_min)
    CALL new_attribute(status, modstr, 'LON_OLD_MIX' &
         , 'valid_max', r = valid_max)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'LEV', p1=LEV, reprid=REPR_LG_CLAMS)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'LEV' &
         , 'creator_of_parameter', c = TRIM(username))
    CALL new_attribute(status, modstr, 'LEV' &
         , 'param_creation_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'LEV' &
         , 'param_modification_time', c = TRIM(cr_date))
    SELECT CASE (TRIM(init_vertcoorname))
    CASE ('press')
       CALL new_attribute(status, modstr, 'LEV' &
            , 'long_name', c = 'Press')
       CALL new_attribute(status, modstr, 'LEV' &
            , 'units', c = 'hPa')
    CASE ('zeta')
       CALL new_attribute(status, modstr, 'LEV' &
            , 'long_name', c = 'Zeta')
       CALL new_attribute(status, modstr, 'LEV' &
            , 'units', c = 'K')
    CASE ('theta')
       CALL new_attribute(status, modstr, 'LEV' &
            , 'long_name', c = 'Theta')
       CALL new_attribute(status, modstr, 'LEV' &
            , 'units', c = 'K')
    CASE DEFAULT
       call error_bi ("Name of vertical coordinate in INIT/POS file not correct!!!",substr)
    END SELECT
    CALL new_attribute(status, modstr, 'LEV' &
         , 'missing_value', r = mdi)
    CALL new_attribute(status, modstr, 'LEV' &
         , 'sample_interval', r = sample_interval)
    CALL new_attribute(status, modstr, 'LEV' &
         , 'flag', c = 'NONE')
    valid_min = 0. 
    valid_max = 100000.
    CALL new_attribute(status, modstr, 'LEV' &
         , 'valid_min', r = valid_min)
    CALL new_attribute(status, modstr, 'LEV' &
         , 'valid_max', r = valid_max)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'LEV_OLD', p1=LEV_OLD, reprid=REPR_LG_CLAMS)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'LEV_OLD' &
         , 'creator_of_parameter', c = TRIM(username))
    CALL new_attribute(status, modstr, 'LEV_OLD' &
         , 'param_creation_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'LEV_OLD' &
         , 'param_modification_time', c = TRIM(cr_date))
    SELECT CASE (TRIM(init_vertcoorname))
    CASE ('press')
       CALL new_attribute(status, modstr, 'LEV_OLD' &
            , 'long_name', c = 'Press')
       CALL new_attribute(status, modstr, 'LEV_OLD' &
            , 'units', c = 'hPa')
    CASE ('zeta')
       CALL new_attribute(status, modstr, 'LEV_OLD' &
            , 'long_name', c = 'Zeta')
       CALL new_attribute(status, modstr, 'LEV_OLD' &
            , 'units', c = 'K')
    CASE ('theta')
       CALL new_attribute(status, modstr, 'LEV_OLD' &
            , 'long_name', c = 'Theta')
       CALL new_attribute(status, modstr, 'LEV_OLD' &
            , 'units', c = 'K')
    CASE DEFAULT
       call error_bi ("Name of vertical coordinate in INIT/POS file not correct!!!",substr)
    END SELECT
    CALL new_attribute(status, modstr, 'LEV_OLD' &
         , 'missing_value', r = mdi)
    CALL new_attribute(status, modstr, 'LEV_OLD' &
         , 'sample_interval', r = sample_interval)
    CALL new_attribute(status, modstr, 'LEV_OLD' &
         , 'flag', c = 'NONE')
    valid_min = 0. 
    valid_max = 100000.
    CALL new_attribute(status, modstr, 'LEV_OLD' &
         , 'valid_min', r = valid_min)
    CALL new_attribute(status, modstr, 'LEV_OLD' &
         , 'valid_max', r = valid_max)
    CALL channel_halt(substr, status)

#ifdef ECHAM5 
! op_sb_20190712+
    ! index position of CLaMS parcels
    CALL new_channel_object(status, modstr, 'POS', p2=pos, &
                reprid=REPR_POS)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'POS' &
         , 'long_name', c='position (index)-lon,lev,lat)' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'POS' &
         , 'units', c=' ' )
    CALL channel_halt(substr, status)
! op_sb_20190712-
#endif

    ! Define channel objects for all parameters 

    DO i = 1, nparams
       if (p_pe==0) &
            write (*,*) 'create channel for ',trim(PARAM(i)%name)
       CALL new_channel_object(status, modstr, trim(PARAM(i)%name), &
            p1=PARAM(i)%values, reprid=REPR_LG_CLAMS)
       CALL channel_halt(substr, status)
    ENDDO

    Do i = 1, nparams_old
       if (p_pe==0) &
            write (*,*) 'create channel for ',trim(PARAM_OLD(i)%name)
       CALL new_channel_object(status, modstr, trim(PARAM_OLD(i)%name), &
            p1=PARAM_OLD(i)%values, reprid=REPR_LG_CLAMS)
       CALL channel_halt(substr, status)
    ENDDO

#ifdef MBM_CLAMS

    write (*,*) 'set variable attributes'
    ! get longname and units for parameters
    DO i = 1, nparams
       PARAM(i)%longname = '???'
       PARAM(i)%units    = '???'
       if (PARAM(i)%name=='EQLAT' .or. PARAM(i)%name=='PV') then
          call nc_get_var_atts (status,pre_thetafile,param(i)%name, &
                             longname=PARAM(i)%longname,units=PARAM(i)%units)
       else
          posloc = INDEX(PARAM(i)%name,'_',back=.true.)
          if (posloc == 0) then
             call nc_get_var_atts (status,pre_metfile,param(i)%name, &
                                longname=PARAM(i)%longname,units=PARAM(i)%units)
          elseif (PARAM(i)%name(posloc+1:)/='TROP1' .and. PARAM(i)%name(posloc+1:)/='TROP2') then
             call nc_get_var_atts (status,pre_metfile,param(i)%name, &
                                longname=PARAM(i)%longname,units=PARAM(i)%units)
          else
              call nc_get_var_atts (status,pre_metfile,param(i)%name(1:posloc-1), &
                                longname=PARAM(i)%longname,units=PARAM(i)%units)
              PARAM(i)%longname = trim(PARAM(i)%longname) // ' at tropopause'
          endif
       endif
    ENDDO

    DO i = 1, nparams_old
       PARAM_OLD(i)%longname = PARAM_OLD(i)%name
       PARAM_OLD(i)%units    = '???'
    ENDDO


#else
    DO i = 1, nparams
       PARAM(i)%longname = PARAM(i)%name
       PARAM(i)%units    = '???'
    ENDDO
    
    DO i = 1, nparams_old
       PARAM_OLD(i)%longname = PARAM_OLD(i)%name
       PARAM_OLD(i)%units    = '???'
    ENDDO
#endif

    ! Set attributes for parameters
    DO i = 1, nparams
       !if (p_pe==0) write (*,*) 'set attributes for ',trim(PARAM(i)%name)
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'creator_of_parameter', c = TRIM(username))
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'param_creation_time', c = TRIM(cr_date))
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'param_modification_time', c = TRIM(cr_date))
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'long_name', c = trim(PARAM(i)%longname) )
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'units', c = trim(PARAM(i)%units))
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'missing_value', r = mdi)
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'sample_interval', r = sample_interval)
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'flag', c = 'NONE')
    ENDDO

    DO i = 1, nparams_old
       !if (p_pe==0) write (*,*) 'set attributes for ',trim(PARAM_OLD(i)%name)
       CALL new_attribute(status, modstr, trim(PARAM_OLD(i)%name) &
            , 'creator_of_parameter', c = TRIM(username))
       CALL new_attribute(status, modstr, trim(PARAM_OLD(i)%name) &
            , 'param_creation_time', c = TRIM(cr_date))
       CALL new_attribute(status, modstr, trim(PARAM_OLD(i)%name) &
            , 'param_modification_time', c = TRIM(cr_date))
       CALL new_attribute(status, modstr, trim(PARAM_OLD(i)%name) &
            , 'long_name', c = trim(PARAM_OLD(i)%longname) )
       CALL new_attribute(status, modstr, trim(PARAM_OLD(i)%name) &
            , 'units', c = trim(PARAM_OLD(i)%units))
       CALL new_attribute(status, modstr, trim(PARAM_OLD(i)%name) &
            , 'missing_value', r = mdi)
       CALL new_attribute(status, modstr, trim(PARAM_OLD(i)%name) &
            , 'sample_interval', r = sample_interval)
       CALL new_attribute(status, modstr, trim(PARAM_OLD(i)%name) &
            , 'flag', c = 'NONE')
    ENDDO


    ! Define channel objects for all clams species
    DO ispec = 1, nspec
       if (p_pe==0 .and. ldiagout) &
            write (*,*) 'create channel for ',trim(SPECARR(ispec)%name)
       CALL new_channel_object(status, modstr, trim(SPECARR(ispec)%name), &
            p1=SPECARR(ispec)%values, reprid=REPR_LG_CLAMS, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, trim(SPECARR(ispec)%name), &
            'ctype',c=trim(SPECARR(ispec)%ctype))
       CALL new_attribute(status, modstr, trim(SPECARR(ispec)%name), &
            'creator_of_parameter', c = trim(username))
       CALL new_attribute(status, modstr, trim(SPECARR(ispec)%name), &
            'param_creation_time',  c = TRIM(cr_date))
       CALL new_attribute(status, modstr, trim(SPECARR(ispec)%name), &
            'param_modification_time', c = TRIM(cr_date))
       CALL new_attribute(status, modstr, trim(SPECARR(ispec)%name), &
            'long_name', c= trim(SPECARR(ispec)%longname))
       CALL new_attribute(status, modstr, trim(SPECARR(ispec)%name), &
            'units', c= trim(SPECARR(ispec)%units))
       CALL new_attribute(status, modstr, trim(SPECARR(ispec)%name), &
            'missing_value', r=mdi )
       CALL new_attribute(status, modstr, trim(SPECARR(ispec)%name), &
            'flag', c='NONE')
       CALL channel_halt(substr, status)
    ENDDO

 ! op_sb_20191011+
#ifdef ECHAM5

   if (nspec+ntrac_cl .gt. maxspec) call error_bi (substr,"To many species: increase MAXSPEC !!!")
   ! get the messy tracer, which should be transported by CLaMS
   DO jt = 1, ntrac_cl
       ispec=nspec+jt
       CALL get_tracer(status,'cl',jt,fullname=name,unit=unit,longname=longname,pxt=ptr_trac)
       SPECARR(ispec)%values => ptr_trac(:,1,1)
       SPECARR(ispec)%name=name
       SPECARR(ispec)%longname=longname
       SPECARR(ispec)%units=unit
       SPECARR(ispec)%ctype=" "
       if (p_pe==0 .and. ldiagout) then
            write (*,*) 'reference channel tracer_cl for ',trim(SPECARR(ispec)%name)
       endif
   enddo
   ! add number of messy tracer (ntrac_cl) to the number of clams tracers'
   nspec=nspec+ntrac_cl
#endif 
  ! op_sb_20191011-

    ! STATE, STATE_VERT (used in MIX and BMIX)
    IF (USE_CLAMSMIX .or. USE_CLAMSBMIX) THEN
       CALL new_channel_object(status, modstr, 'STATE', p1=STATE, reprid=REPR_LG_CLAMS)
       CALL channel_halt(substr, status)
       CALL new_channel_object(status, modstr, 'STATE_VERT', p1=STATE_VERT, &
            reprid=REPR_LG_CLAMS)
       CALL channel_halt(substr, status)
    ENDIF


    ! THETA_OLD_MIX and BVF_OLD_MIX (for vertical mixing only)
    IF (USE_CLAMSMIX .and. switch_mixing==2) THEN 
       CALL new_channel_object(status, modstr, 'THETA_OLD_MIX', &
            p1=THETA_OLD_MIX, reprid=REPR_LG_CLAMS)
       CALL channel_halt(substr, status)
       CALL new_channel_object(status, modstr, 'BVF_OLD_MIX', &
            p1=BVF_OLD_MIX, reprid=REPR_LG_CLAMS)
       CALL channel_halt(substr, status)
    ENDIF

    ! LEV_GRID, LEV_DELTA and R_GRID (used in MIX and BMIX)
    IF (USE_CLAMSMIX .or. USE_CLAMSBMIX) THEN
       CALL new_channel_object(status, modstr, 'LEV_GRID', p1=LEV_GRID, &
            reprid=REPR_MIX_GRID)
       CALL channel_halt(substr, status)
       CALL new_channel_object(status, modstr, 'LEV_DELTA', p1=LEV_DELTA, &
            reprid=REPR_MIX_GRID)
       CALL channel_halt(substr, status)
       IF (ASSOCIATED(adapt_par%r_grid)) THEN
          CALL new_channel_object(status, modstr, 'R_GRID', p1=R_GRID, &
               reprid=REPR_MIX_GRID)
          CALL channel_halt(substr, status)
       ENDIF
    ENDIF

  
    !-----------------------------------------------------------------
    ! Channel CLAMS: add ECHAM5 channel objects
    !-----------------------------------------------------------------

#ifdef ECHAM5
    CALL new_channel_object(status, modstr, 'E5_ZETA', p3=E5_ZETA_D, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'E5_ZETA' &
         , 'creator_of_parameter', c = TRIM(username))
    CALL new_attribute(status, modstr, 'E5_ZETA' &
         , 'param_creation_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'E5_ZETA' &
         , 'param_modification_time', c = TRIM(cr_date))
    CALL new_attribute(status, modstr, 'E5_ZETA' &
         , 'long_name', c = TRIM(init_vertcoorname))
    SELECT CASE (TRIM(init_vertcoorname))
       CASE ('press')
          valid_min = 0. 
          valid_max = 100000.
       CASE ('zeta')
          valid_min = 0.
          valid_max = 8000.
       CASE ('theta')
          valid_min = 0. 
          valid_max = 3000.
       CASE DEFAULT
          call error_bi ("Name of vertical coordinate in INIT/POS file not correct!!!",substr)
       END SELECT
    CALL new_attribute(status, modstr, 'E5_ZETA' &
         , 'valid_min', r = valid_min)
    CALL new_attribute(status, modstr, 'E5_ZETA' &
         , 'valid_max', r = valid_max)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'E5_LON3D', p3=E5_LON3D_D, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr, 'E5_LAT3D', p3=E5_LAT3D_D, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr, 'E5_OMEGA', p3=E5_OMEGA_D, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr, 'E5_THETADOT', p3=E5_THETADOT_D, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr, 'E5_THETADOT2', p3=E5_THETADOT2_D, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr, 'E5_ZETADOT', p3=E5_ZETADOT_D, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'E5_SIGMADOT', p3=E5_SIGMADOT_D, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr, 'E5_SIGMA', p3=E5_SIGMA_D, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'E5_PSDOT', p2=E5_PSDOT_D, &
         reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr, 'E5_HEATING', p3=E5_HEATING_D, reprid=GP_3D_MID)
    CALL new_attribute(status, modstr, 'E5_HEATING' &
         , 'long_name', c = 'diabatic heating rate Q')
    CALL channel_halt(substr, status)
    ! op_pj_20180614+
    CALL new_attribute(status, modstr, 'E5_HEATING' &
         , 'units', c = 'K/s')
    CALL channel_halt(substr, status)
    ! op_pj_20180614-
! op_sb_20190730+
    CALL new_channel_object(status, modstr, 'SPC_G' &
         , p3=spc_g, reprid=GP_3D_MID, lrestreq=.FALSE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'SPC_G' &
         , 'long_name', c='number of cells in each grid box' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'SPC_G' &
         , 'units', c=' ' )
    CALL channel_halt(substr, status)
! op_sb_20190730-
#endif


    !-----------------------------------------------------------------
    ! Define channel winddata
    !-----------------------------------------------------------------

    CALL new_channel(status, 'winddata', lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'winddata', 'winddata'//'_version', c = modver)
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, 'winddata', 'UDT', p3=UDT, reprid=REPR_3DINP_CLAMS)
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, 'winddata', 'VDT', p3=VDT, reprid=REPR_3DINP_CLAMS)
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, 'winddata', 'WDT', p3=WDT, reprid=REPR_3DINP_CLAMS)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, 'winddata', 'LEVELDT', p3=leveldt, &
         reprid=REPR_3DINP_CLAMS)
    CALL channel_halt(substr, status)
 
    CALL new_channel_object(status, 'winddata', 'DLEVDZDT', p3=dlevdzdt, &
         reprid=REPR_3DINP_CLAMS)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, 'winddata', 'UFUT', p3=UFUT, reprid=REPR_3DINP_CLAMS)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, 'winddata', 'VFUT', p3=VFUT, reprid=REPR_3DINP_CLAMS)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, 'winddata', 'WFUT', p3=WFUT, reprid=REPR_3DINP_CLAMS)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, 'winddata', 'LEVELFUT', p3=levelfut, &
         reprid=REPR_3DINP_CLAMS)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, 'winddata', 'DLEVDZFUT', p3=dlevdzfut, &
         reprid=REPR_3DINP_CLAMS)
    CALL channel_halt(substr, status)

    IF (USE_CLAMSSEDI) THEN

       CALL new_channel_object(status, 'winddata', 'UDT_sedi', p3=UDT_sedi, &
            reprid=REPR_3DINP_CLAMS)
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, 'winddata', 'VDT_sedi', p3=VDT_sedi, &
            reprid=REPR_3DINP_CLAMS)
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, 'winddata', 'WDT_sedi', p3=WDT_sedi, &
            reprid=REPR_3DINP_CLAMS)
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, 'winddata', 'LEVELDT_sedi', p3=leveldt_sedi, &
            reprid=REPR_3DINP_CLAMS)
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, 'winddata', 'DLEVDZDT_sedi', p3=dlevdzdt_sedi, &
            reprid=REPR_3DINP_CLAMS)
       CALL channel_halt(substr, status)


    ENDIF


    DO i = 1, nparams

#ifdef MBM_CLAMS       
       if (PREDATA(i)%name=='EQLAT' .or. PREDATA(i)%name=='PV') then

          if (p_pe==0 .and. ldiagout) &
               write (*,*) 'create channel for ',"PREDATA_",trim(PREDATA(i)%name)
          CALL new_channel_object(status, 'winddata', "PREDATA_"//trim(PREDATA(i)&
               &%name), p3=PREDATA(i)%values, reprid=REPR_3DINP_CLAMSTHETA)
          CALL channel_halt(substr, status)
          
          if (p_pe==0 .and. ldiagout) &
               write (*,*) 'create channel for ',"FUTDATA_",trim(PREDATA(i)%name)
          CALL new_channel_object(status, 'winddata', "FUTDATA_"//trim(FUTDATA(i)&
               &%name), p3=FUTDATA(i)%values, reprid=REPR_3DINP_CLAMSTHETA)
          CALL channel_halt(substr, status)
#else
       if (PREDATA(i)%name=='EQLAT') then
             !!! ???
       elseif (PREDATA(i)%name=='PV') then

          if (p_pe==0 .and. ldiagout) &
               write (*,*) 'create channel for ',"PREDATA_",trim(PREDATA(i)%name)
          CALL new_channel_object(status, 'winddata', "PREDATA_"//trim(PREDATA(i)&
               &%name), p3=PREDATA(i)%values, reprid=REPR_3DINP_CLAMS)
          CALL channel_halt(substr, status)
          
          if (p_pe==0 .and. ldiagout) &
               write (*,*) 'create channel for ',"FUTDATA_",trim(PREDATA(i)%name)
          CALL new_channel_object(status, 'winddata', "FUTDATA_"//trim(FUTDATA(i)&
               &%name), p3=FUTDATA(i)%values, reprid=REPR_3DINP_CLAMS)
          CALL channel_halt(substr, status)
         
#endif          
          
       else

          if (p_pe==0 .and. ldiagout) &
               write (*,*) 'create channel for ',"PREDATA_",trim(PREDATA(i)%name)
          CALL new_channel_object(status, 'winddata', "PREDATA_"//trim(PREDATA(i)&
               &%name), p3=PREDATA(i)%values, reprid=REPR_3DINP_CLAMS)
          CALL channel_halt(substr, status)
          
          if (p_pe==0 .and. ldiagout) &
               write (*,*) 'create channel for ',"FUTDATA_",trim(PREDATA(i)%name)
          CALL new_channel_object(status, 'winddata', "FUTDATA_"//trim(FUTDATA(i)&
               &%name), p3=FUTDATA(i)%values, reprid=REPR_3DINP_CLAMS)
          CALL channel_halt(substr, status)

       endif
       
       PREDATA(i)%longname = '???'
       PREDATA(i)%units    = '???'
#ifdef MBM_CLAMS
       if (PREDATA(i)%name=='EQLAT' .or. PREDATA(i)%name=='PV') then
          call nc_get_var_atts (status,pre_thetafile,PREDATA(i)%name, &
                             longname=PREDATA(i)%longname,units=PREDATA(i)%units)
      else
         posloc = INDEX(PARAM(i)%name,'_',back=.true.)
         if (posloc == 0) then
            call nc_get_var_atts (status,pre_metfile,PREDATA(i)%name, &
                                  longname=PREDATA(i)%longname,units=PREDATA(i)%units)
          elseif (PARAM(i)%name(posloc+1:)/='TROP1' .and. PARAM(i)%name(posloc+1:)/='TROP2') then
            call nc_get_var_atts (status,pre_metfile,PREDATA(i)%name, &
                                  longname=PREDATA(i)%longname,units=PREDATA(i)%units)
          else
             call nc_get_var_atts (status,pre_metfile,PREDATA(i)%name(1:posloc-1), &
                                  longname=PREDATA(i)%longname,units=PREDATA(i)%units)
             PREDATA(i)%longname = trim(PREDATA(i)%longname) // ' at tropopause'
          endif
       endif
#endif
       FUTDATA(i)%longname = PREDATA(i)%longname
       FUTDATA(i)%units    = PREDATA(i)%units

    ENDDO

#ifdef ECHAM5

! ju_ec_20180717+
! Other coupling for CLaMS parcel gridding
! Set ECHAM grid size variables

  ! Allocate position array --> now a channel object
!!$  ALLOCATE(pos(3,dnparts_max))   !  op_sb_20190712

  ! ju_ec_20180717 I needed the E5_LON/LAT/LEVEL lines so that I could obtain
  ! the EMAC grid shape in preparation for EMAC-CLaMS coupling. These values
  ! are also obtained in clams_global_end, and are re-acquired on each call. I 
  ! don't see why, as they shouldn't change, should they? Perhaps 
  ! clams_global_end could do for some more cleaning.
  CALL get_channel_object_dimvar(status, 'g2a', 'um1', dva, units)
  E5_LON   => dva(1)%ptr(:)
  E5_LAT   => dva(3)%ptr(:)
  E5_LEVEL => dva(2)%ptr(:)

  nlat  = size(E5_LAT  )
  nlon  = size(E5_LON  )
  nlev  = size(E5_LEVEL)
  nlev_i = nlev ! op_sb_20191028

  ! Allocate boundary information
  ALLOCATE(E5_LAT_BOUNDS(nlat+1))

  ALLOCATE(E5_LON_BOUNDS(nlon  ))

  ! Allocate parcel count arrays (parcels per grid box)
  ALLOCATE( pc_d(nlon,nlev,nlat) )  ! local/decomposed
  ALLOCATE( pc_g(nlon,nlev,nlat) )  ! global
! ju_ec_20180717-
  
  if (clams_gridding) then
      ALLOCATE( clgp_data(nlon,nlev,n_cltr,nlat) )
      clgp_data(:,:,:,:) = 0.0_dp
  endif
#endif
  
 END SUBROUTINE clams_init_memory
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clams_init_coupling(flag)

#ifdef ECHAM5
    USE messy_main_blather_bi,   ONLY: start_message_bi, end_message_bi
    USE messy_main_mpi_bi,       ONLY: p_pe
    ! op_sb_20191022+
    USE messy_main_transform_bi, ONLY: trp_gpdc_gpgl
#endif
    USE messy_main_blather_bi,   ONLY: error_bi
    USE messy_main_channel,      ONLY: get_channel_object, new_attribute, get_attribute, &
                                       new_channel_object_reference ! op_sb_20191011
    USE messy_main_channel_error_bi,  ONLY: channel_halt

    USE messy_clams_global,      ONLY: username, mdi, nparams, sample_interval, &
                                       PARAM, PREDATA, FUTDATA, &
                                       E5_LAT, E5_LON, E5_LEVEL, & !EMAC-CLaMS Coupl.
                                       n_cltr, cl_grid_tracers, &
                                       clams_gridding, clams_grid_verbose, rank, &
                                       nspec, SPECARR, ldiagout! op_sb_20191014
    ! op_sb_20191018+
!    use messy_clamsmix_global,   ONLY: nmixspec
    ! op_sb_20191018-

    
! op_pj_20180614+
#ifdef ECHAM5
#ifdef MESSYTENDENCY
    USE messy_main_tendency_bi,  ONLY: mtend_request, mtend_id_t
    USE messy_main_blather_bi,   ONLY: error_bi
#endif
#endif
! op_pj_20180614-

#ifdef ECHAM5
! ju_ec_20180627+
! USE commands for CLaMS tracers
    USE messy_main_tracer_mem_bi,   ONLY: CLTRSTR, &
                                          ntrac_cl  ! op_sb_20191014
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_main_tracer,          ONLY: get_tracer
! ju_ec_20180627-
#endif
    
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: flag
    

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr =  'clams_init_coupling'
    INTEGER :: status, i, jt
    iNTEGER :: iLat, iLon ! ju_ec_20180717 for EMAC boundary determination
    INTEGER :: ispec ! op_sb_20191014

    CHARACTER(40) :: cr_date
    CHARACTER(8)  :: ydate
    CHARACTER(10) :: ytime
    CHARACTER(8)  :: name ! op_sb_20191018

    SELECT CASE(flag)

    CASE (1)
#ifdef ECHAM5
    CALL start_message_bi(modstr, 'COUPLING TO DRIVER FIELDS', substr)

    CALL get_channel_object(status, 'g2a', 'um1', p3=E5_UWIND_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'g2a', 'vm1', p3=E5_VWIND_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'ECHAM5', 'tm1', p3=E5_TEMP_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'ECHAM5', 'press', p3=E5_PRESS_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'ECHAM5', 'tpot', p3=E5_THETA_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'g3b', 'aps', p2=E5_PSURF_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'scnbuf', 'vervel', p3=E5_OMEGA_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'scnbuf', 'tte', p3=E5_TTE_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'scnbuf', 'alpste', p2=E5_ALPSTE_D)
    CALL channel_halt(substr, status)
! op_pj_20180614+
!!$    CALL get_channel_object(status, 'tdiag', 'dtdt_vdiff', p3=E5_HEAT_VDIFF_D)
!!$    CALL channel_halt(substr, status)
!!$    CALL get_channel_object(status, 'tdiag', 'dtdt_cloudconv', p3=E5_HEAT_CLOUDCONV_D)
!!$    CALL channel_halt(substr, status)
!!$    CALL get_channel_object(status, 'tdiag', 'dtdt_rheat', p3=E5_HEAT_RHEAT_D)
!!$    CALL channel_halt(substr, status)
!!$    CALL get_channel_object(status, 'tdiag', 'dtdt_sso', p3=E5_HEAT_SSO_D)
!!$    CALL channel_halt(substr, status)
! op_pj_20180614-
    CALL get_channel_object(status, 'ECHAM5', 'qm1', p3=E5_SH_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'ECHAM5', 'xim1', p3=E5_IWC_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'ECHAM5', 'xlm1', p3=E5_CLWC_D)
    CALL channel_halt(substr, status)
    ! op_sb_20190822+
    CALL get_channel_object(status, 'tropop', 'pblh_i', p2=pblh_i)
    CALL channel_halt(substr, status)

    ! dry air mass for transformation routines
    CALL get_channel_object(status, 'gridd_def', 'grmassdry', p3=grmass)
    CALL channel_halt(substr, status)
  
    ! point from tracer_cl memory to SPECARR clams tracer (avoid to allocate memory)
    DO jt = nspec - ntrac_cl+1, nspec
       CALL new_channel_object_reference(status, 'tracer_cl',&
                                  trim(SPECARR(jt)%name),modstr,trim(SPECARR(jt)%name))
       CALL channel_halt(substr, status)
    enddo
    ! op_sb_20191022-
#endif

      
#ifndef MBM_CLAMS

    ! get longname and units for parameters
    DO i = 1, nparams
       SELECT CASE (PARAM(i)%name)
       CASE ('TEMP')
          CALL get_attribute (status,'ECHAM5', 'tm1', 'long_name',c=PARAM(i)%longname)
          if (p_pe==0) write (*,*) 'ECHAM_TEMP: long_name=',trim (PARAM(i)%longname)
          CALL get_attribute (status,'ECHAM5', 'tm1', 'units',c=PARAM(i)%units)
          if (p_pe==0) write (*,*) 'ECHAM_TEMP: units=',trim(PARAM(i)%units)
       CASE ('PRESS')
          CALL get_attribute (status,'ECHAM5', 'press', 'long_name',c=PARAM(i)%longname)
          if (p_pe==0) write (*,*) 'ECHAM_TEMP: long_name=',trim (PARAM(i)%longname)
          CALL get_attribute (status,'ECHAM5', 'press', 'units',c=PARAM(i)%units)
          if (p_pe==0) write (*,*) 'ECHAM_TEMP: units=',trim(PARAM(i)%units)
          ! Change units of pressure (Pa -> hPa)
          PARAM(i)%units    = 'hPa'
       CASE ('H2O')
          CALL get_attribute (status,'ECHAM5', 'qm1', 'long_name',c=PARAM(i)%longname)
          if (p_pe==0) write (*,*) 'ECHAM_TEMP: long_name=',trim (PARAM(i)%longname)
          CALL get_attribute (status,'ECHAM5', 'qm1', 'units',c=PARAM(i)%units)
          if (p_pe==0) write (*,*) 'ECHAM_TEMP: units=',trim(PARAM(i)%units)
          ! Conleved from mass mixing ratio to volume mixing ratio
          PARAM(i)%longname = 'H2O volume mixing ratio'
          PARAM(i)%units    = 'm**3/m**3'
       CASE ('IWC')
          CALL get_attribute (status,'ECHAM5', 'xim1', 'long_name',c=PARAM(i)%longname)
          if (p_pe==0) write (*,*) 'ECHAM_TEMP: long_name=',trim (PARAM(i)%longname)
          CALL get_attribute (status,'ECHAM5', 'xim1', 'units',c=PARAM(i)%units)
          if (p_pe==0) write (*,*) 'ECHAM_TEMP: units=',trim(PARAM(i)%units)
       CASE ('CLWC')
          CALL get_attribute (status,'ECHAM5', 'xlm1', 'long_name',c=PARAM(i)%longname)
          if (p_pe==0) write (*,*) 'ECHAM_TEMP: long_name=',trim (PARAM(i)%longname)
          CALL get_attribute (status,'ECHAM5', 'xlm1', 'units',c=PARAM(i)%units)
          if (p_pe==0) write (*,*) 'ECHAM_TEMP: units=',trim(PARAM(i)%units)
       CASE ('ZETA')
          PARAM(i)%longname = 'Hybrid Theta'
          PARAM(i)%units    = 'K'
       CASE ('THETA')
          CALL get_attribute (status,'ECHAM5', 'tpot', 'long_name',c=PARAM(i)%longname)
          if (p_pe==0) write (*,*) 'ECHAM_TEMP: long_name=',trim (PARAM(i)%longname)
          CALL get_attribute (status,'ECHAM5', 'tpot', 'units',c=PARAM(i)%units)
          if (p_pe==0) write (*,*) 'ECHAM_TEMP: units=',trim(PARAM(i)%units)
       CASE DEFAULT
          PARAM(i)%longname = '???'
          PARAM(i)%units    = '???'
       END SELECT

       PREDATA(i)%longname = PARAM(i)%longname
       PREDATA(i)%units    = PARAM(i)%units
       FUTDATA(i)%longname = PARAM(i)%longname
       FUTDATA(i)%units    = PARAM(i)%units

   ENDDO

    CALL DATE_AND_TIME(ydate, ytime)
    WRITE(cr_date,'(A4,5(A,A2))')   &
         ydate(1:4),'-',ydate(5:6),'-',ydate(7:8),' ',  &
         ytime(1:2),':',ytime(3:4),':',ytime(5:6) 

    ! Set attributes for parameters
    DO i = 1, nparams
       !if (p_pe==0) write (*,*) 'set attributes for ',trim(PARAM(i)%name)
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'creator_of_parameter', c = TRIM(username))
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'param_creation_time', c = TRIM(cr_date))
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'param_modification_time', c = TRIM(cr_date))
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'long_name', c = trim(PARAM(i)%longname) )
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'units', c = trim(PARAM(i)%units))
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'missing_value', r = mdi)
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'sample_interval', r = sample_interval)
       CALL new_attribute(status, modstr, trim(PARAM(i)%name) &
            , 'flag', c = 'NONE')

    ENDDO
#endif

#ifdef ECHAM5

! op_pj_20180614+
!!$#ifdef MESSYTENDENCY
!!$       CALL mtend_request('dyn', mtend_id_t, tte_dyn)
!!$#endif
! op_pj_20180614-

    CALL end_message_bi(modstr, 'COUPLING TO DRIVER FIELDS', substr)
#endif
     
   

#ifdef ECHAM5

! ju_ec_20190206+
    if ((clams_gridding).AND.(n_cltr.EQ.0)) &
        CALL error_bi( &
                  "CLaMS gridding is on, but there are no tracers!"&
                , substr&
                )
! ju_ec_20190206-

!
! ju_ec_20190206
!
! Commented to move to clams_new_tracer. This makes more sense there, because
! the tracer indexes shouldn't change. Also, I need to know what the indexes are
! in that routine, so I can assign properties to them via SET_TRACER from the
! messy main tracer routines.
!
!! ju_ec_20180627+
!! Coupling to CLaMS tracers
!    if (clams_gridding) then
!        if (n_cltr .EQ. 0) write(*,*) &
!            "CLaMS gridding is on, but there are no tracers!"
!        CLTR_IDX(:) = 0
!        DO i = 1, n_cltr
!           CALL get_tracer( status, CLTRSTR, cl_grid_tracers(i)%name, idx=CLTR_IDX(i) )
!           IF ((rank.eq.0).and.(clams_grid_verbose)) write(*,*)&
!               "Gave tracer ",cl_grid_tracers(i)%name," index ",CLTR_IDX(i)
!           IF (status /= 0) CLTR_IDX(i) = 0
!        END DO
!
!
!    endif
!! ju_ec_20180627-
  
! ju_ec_20180717+
! Other coupling for CLaMS parcel gridding
! Set ECHAM grid size variables

  ! Determine boundary latitudes
  E5_LAT_BOUNDS(1) = 90.
  DO iLat = 1, nlat-1
    E5_LAT_BOUNDS(iLat+1) = (E5_LAT(iLat+1)+E5_LAT(iLat))/2.
  END DO
  E5_LAT_BOUNDS(nlat+1) = -90.

  ! Determine boundary longitudes
  DO iLon = 1, nlon-1
    E5_LON_BOUNDS(iLon) = (E5_LON(iLon+1)+E5_LON(iLon))/2.
  END DO
   E5_LON_BOUNDS(nlon) = (E5_LON(nlon)+360.)/2.
! ju_ec_20180717-

 
#endif

   CASE (2)
#ifdef ECHAM5
#ifdef MESSYTENDENCY
       CALL mtend_request('dyn', mtend_id_t, tte_dyn)
#endif
#endif
    CASE DEFAULT
       CALL error_bi('UNKNOWN FLAG FOR CALL ',substr)
    END SELECT
       
  END SUBROUTINE clams_init_coupling
!-----------------------------------------------------------------------

! op_pj_20180614+
!---------------------------------------------------------------------------
  SUBROUTINE clams_local_start

#ifdef ECHAM5
    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi, ONLY: jrow
    USE messy_main_data_bi,         ONLY: tte_3d

    IMPLICIT NONE

#ifndef MESSYTENDENCY
    ! save temperature tendency after advection
    E5_HEATING_D(:,:,jrow) = tte_3d(:,:,jrow)
#endif

#endif

  END SUBROUTINE clams_local_start
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE clams_local_end

#ifdef ECHAM5
    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi, ONLY: jrow
    USE messy_main_data_bi,         ONLY: tte_3d
    
    IMPLICIT NONE

    ! calculate dtdt for diabatic processes only
#ifndef MESSYTENDENCY
    E5_HEATING_D(:,:,jrow) = tte_3d(:,:,jrow) - E5_HEATING_D(:,:,jrow)
#else
    E5_HEATING_D(:,:,jrow) = tte_3d(:,:,jrow) - tte_dyn(:,:,jrow)
#endif

#endif

  END SUBROUTINE clams_local_end
!---------------------------------------------------------------------------
! op_pj_20180614-

!-----------------------------------------------------------------------
  SUBROUTINE clams_global_start
    ! Read INIT file and meteorological data
    ! collect ECHAM wind fields if ECHAM5 is base model

    ! BMIL
    USE messy_main_mpi_bi,      ONLY: p_pe, p_sum
    USE messy_main_blather_bi,  ONLY: error_bi
    USE messy_main_channel_error_bi,  ONLY: channel_halt
    USE messy_main_channel_bi,  ONLY: main_channel_write_output
    USE messy_main_timer_bi,    ONLY: event_state

    ! SMIL: CLaMS
    USE messy_clamsmix_si,     ONLY: mixevent
    USE messy_clamsbmix_si,    ONLY: bmixevent
    USE messy_clamscirrus_si,  ONLY: cirrusevent
!!#D clamschem +
    USE messy_clamschem_si,    ONLY: chemevent
!!#D clamschem -
    USE messy_clamssedi_si,    ONLY: sedievent
    USE messy_clamstracer_si,  ONLY: tracerevent
    USE messy_clamsdeepconv_si,ONLY: deepconvevent

    ! SMCL: MESSy
    USE messy_main_timer,      ONLY: lstart, lresume, current_date,          &
                                     YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, &
                                     YEAR_START, MONTH_START, DAY_START,     &
                                     HOUR_START, MINUTE_START, SECOND_START, &
                                     YEAR_NEXT, MONTH_NEXT, DAY_NEXT,        &
                                     HOUR_NEXT, MINUTE_NEXT, SECOND_NEXT,    &
                                     delta_time, lfirst_cycle
    USE messy_main_switch,     ONLY: USE_CLAMSTRAJ, USE_CLAMSMIX, &
                                     USE_CLAMSBMIX,   &
                                     USE_CLAMSCIRRUS, &
                                     USE_DISSOC, USE_CLAMS, &
                                     USE_CLAMSSEDI, USE_CLAMSTRACER, &
                                     USE_CLAMSDEEPCONV
!!#D clamschem +
    USE messy_main_switch,     ONLY: USE_CLAMSCHEM, USE_CLAMSCHEME5
!!#D clamschem -
    USE messy_main_channel,    ONLY: set_channel_output
    ! IF  clamssedi would be a correct sub-submodel of clams this 
    ! MESSy non-conformuse use would not be necessary.
    USE messy_main_switch,     ONLY: USE_CLAMSSEDI

    ! SMCL: CLaMS
    USE messy_clams_global,    ONLY: DP, prec, mdi, eps, &
                                     initfile, init_vertcoorname, &
                                     nspec, SPECARR,  &
                                     dnparts, dnparts_max, dnparts_co, dnparts_max_co, &
                                     nparts, idx, &
                                     lmixevent, lbmixevent, lcirrusevent, &
                                     lchemevent, lsedievent, ltracerevent, &
                                     ldeepconvevent, lclamsoutevent, &
                                     init_h2o_emac,&
                                     lperpetuum, ldiagout, &
                                     pre_year, pre_month, pre_day, pre_sec, &
                                     pre_year_co, pre_month_co, pre_day_co, pre_sec_co, &
                                     irdday, irdsec, &
                                     fut_year, fut_month, fut_day, fut_sec, &
                                     buffersize
    USE messy_clams_global,    ONLY: CYR=> YEAR, CMM=> MONTH, CDY=>DAY      &
                                   , CHR=> HOUR, CMI=>MINUTE, CSE=>SECOND   &
                                   , CYRN=>YEAR_NEXT,    CMMN=>MONTH_NEXT   &
                                   , CDYN=>DAY_NEXT,     CHRN=>HOUR_NEXT    &
                                   , CMIN=>MINUTE_NEXT,  CSEN=>SECOND_NEXT  &
                                   , CYRs=>YEAR_START,   CMMs=>MONTH_START  &
                                   , CDYs=>DAY_START,    CHRs=>HOUR_START   &
                                   , CMIs=>MINUTE_START, CSEs=>SECOND_START &
                                   , cdelta_time => delta_time              &
                                   , clstart => lstart, clres => lresume    &
                                   , cl1cyc => lfirst_cycle

!!#D clamschem +
    USE messy_clamschem_global,    ONLY: timestep_chem, rates, const
!!#D clamschem -
    USE messy_clamsbmix_global,    ONLY: timestep_bmix
    USE messy_clamsmix_global,     ONLY: timestep_mix
    USE messy_clamssedi_global,    ONLY: timestep_sedi
    USE messy_clamstracer_global,  ONLY: timestep_tracer
    USE messy_clamsdeepconv_global,ONLY: timestep_deepconv
    USE messy_clams_read_metdata,  ONLY: clams_read_ecmwf_files
    USE messy_clams_tools_utils,   ONLY: uppercase
    USE messy_dissoc,              ONLY: timestep_dissoc

    USE messy_main_tracer_mem_bi,  ONLY: ntrac_cl
    USE netcdf

    IMPLICIT NONE

    REAL(PREC), ALLOCATABLE, DIMENSION(:) :: field_init
    REAL(DP),DIMENSION(1)                 :: help_array
    INTEGER :: ncid, rcode, varid
    INTEGER :: i, ipart
    INTEGER :: status

    CHARACTER(LEN=*), PARAMETER :: substr = 'clams_global_start'
    LOGICAL :: init_sh
    
    init_sh = .FALSE.

    IF (p_pe==0 .and. ldiagout) write(*,*) uppercase(substr)

    if (lresume) then
       dnparts = dnparts_co(1)
       dnparts_max = dnparts_max_co(1)
       nparts  = p_sum(dnparts)
       write(*,*) 'nach restart: p_pe', p_pe, ' dnparts ', dnparts,' nparts',nparts
    else
       dnparts_co(1) = dnparts
       dnparts_max_co(1) = dnparts_max
    endif

#ifdef MBM_CLAMS   
    if(lstart) then
       pre_year  = YEAR_START
       pre_month = MONTH_START
       pre_day   = DAY_START
       pre_sec   = HOUR_START*3600 + MINUTE_START*60 + SECOND_START

       pre_year_co  = YEAR_START
       pre_month_co = MONTH_START
       pre_day_co   = DAY_START
       pre_sec_co   = HOUR_START*3600 + MINUTE_START*60 + SECOND_START
    elseif (lresume) then
       pre_year  = pre_year_co  
       pre_month = pre_month_co 
       pre_day   = pre_day_co   
       pre_sec   = pre_sec_co  
    end if
#endif

    
    !----------------------------------------------------------------------------
    !  check timesteps of submodels:
    !
    !  - timestep_mix must be equal to or multiple of timestep_chem
    !     (LAT_OLD etc are updated in mix) !
    !      
    !  - timestep_bmix must be equal to or multiple of timestep_chem
    !     (LAT_OLD etc are updated in bmix) !
    !      
    !  - timestep_dissoc must be equal to timestep_chem
    !     -> use correct photolysis rates
    !
    !----------------------------------------------------------------------------
    if (lstart) then
       
       if (USE_CLAMSCHEM) then
          
         if (USE_CLAMSMIX) then
             if (mod(timestep_mix*3600,timestep_chem) /= 0) then
                call error_bi &
                     ("Wrong timesteps: mixing timestep must be equal to or multiple of chemistry timestep !!!",substr)
             endif
          endif
          if (USE_CLAMSBMIX) then
             if (mod(timestep_bmix*3600,timestep_chem) /= 0) then
                call error_bi &
                     ("Wrong timesteps: bmix timestep must be equal to or multiple of chemistry timestep !!!",substr)
             endif
          endif
          
          if (USE_DISSOC) then
             if (timestep_dissoc /= timestep_chem) then
                call error_bi &
                     ("Wrong timesteps: dissoc timestep must be equal to chemistry timestep !!!",substr)
             endif
          endif

          
       endif
       
    endif       ! timestep_mix/timestep_bmix must be equal to or multiple of timestep_chem
       ! (LAT_OLD etc are updated in mix/bmix) !



    ! set date for internal usage in core
    CYR = YEAR
    CMM = MONTH
    CDY = DAY
    CHR = HOUR
    CMI = MINUTE
    CSE = SECOND
    CYRN = YEAR_NEXT
    CMMN = MONTH_NEXT
    CDYN = DAY_NEXT
    CHRN = HOUR_NEXT
    CMIN = MINUTE_NEXT
    CSEN = SECOND_NEXT
    CYRS = YEAR_START
    CMMS = MONTH_START
    CDYS = DAY_START
    CHRS = HOUR_START
    CMIS = MINUTE_START
    CSES = SECOND_START
    cdelta_time = delta_time
    clstart = lstart
    cl1cyc  = lfirst_cycle
    clres   = lresume


    !*******************************************************
    ! Read meteorological dataset
    !*******************************************************
#ifdef MBM_CLAMS
    !if (USE_CLAMSTRAJ .or. USE_CLAMSSEDI) then
       !if (p_pe==0) write (*,*) 'read ecmwf file'
       CALL clams_read_ecmwf_files (status, USE_CLAMSSEDI)
       if (status /= 0 ) call error_bi ("Cannot read ECMWF windfile",substr)
    !endif

       ! in clams_read_ecmwf_files the following variables are assigned:
       ! - udt, vdt, wdt, leveldt
       ! - ufut, vfut, wfut, levelfut
       ! - PREDATA, FUTDATA

       
#endif
 
#ifdef ECHAM5
    !--------------------------------------------------------------
    ! Interval between met data
    !--------------------------------------------------------------
    ! CLAMS: in messy_clams_read_ecmwf 
    irdday = 0
    irdsec = delta_time

    ! The following variables are assigned in clams_global_end:
    ! - udt, vdt, wdt, leveldt
    ! - ufut, vfut, wfut, levelfut
    ! - PREDATA, FUTDATA

#endif

    !*******************************************************
    ! Read INIT-File only in the first cycle:
    !*******************************************************
    if (lstart) then

       IF (p_pe==0) write(*,*) uppercase(substr), ': in lstart'
 
       ! Open init file
       rcode = nf90_open(initfile, nf90_nowrite, ncid, buffersize)
       if (rcode /= nf90_noerr) call error_bi(substr,"Cannot open file "//trim(initfile))

       !--------------------------------------------------------------
       ! Read initial positions from init file
       !--------------------------------------------------------------

       if (p_pe==0) WRITE (*,*)
       if (p_pe==0) WRITE (*,*) 'Read initial positions from ', trim(initfile)
       if (p_pe==0) WRITE (*,*)
     
       ! read time
       !write (*,*) 'read time'
       rcode = nf90_inq_varid (ncid,'time',varid)
       IF (rcode /= 0) CALL error_bi('Cannot find variable time ',substr)
       rcode = nf90_get_var (ncid,varid,help_array)
       IF (rcode /= 0) CALL error_bi('Cannot read variable time ',substr)
       JULTIME = help_array(1)

       ! read latitudes
       !write (*,*) 'read lat'
       rcode = nf90_inq_varid (ncid,'LAT',varid)
       IF (rcode /= 0) CALL error_bi('Cannot find variable LAT ',substr)
       rcode = nf90_get_var (ncid,varid,LAT(1:dnparts),start=(/idx(p_pe,1)/), count=(/dnparts/))
       IF (rcode /= 0) CALL error_bi('Cannot read variable LAT ',substr)
       LAT(dnparts+1:dnparts_max) = mdi

       ! read longitude
       !write (*,*) 'read lon'
       rcode = nf90_inq_varid (ncid,'LON',varid)
       IF (rcode /= 0) CALL error_bi('Cannot find variable LON ',substr)
       rcode = nf90_get_var (ncid,varid,LON(1:dnparts),start=(/idx(p_pe,1)/), count=(/dnparts/))
       IF (rcode /= 0) CALL error_bi('Cannot read variable LON ',substr)
       LON(dnparts+1:dnparts_max) = mdi 

       ! read level
#ifdef MBM_CLAMS
       !write (*,*) 'read ',TRIM(uppercase(init_vertcoorname))
       rcode = nf90_inq_varid (ncid,TRIM(uppercase(init_vertcoorname)),varid)
       IF (rcode /= 0) &
            CALL error_bi('Cannot find variable '//TRIM(uppercase(init_vertcoorname)),substr)
       rcode = nf90_get_var(ncid,varid,LEV(1:dnparts),start=(/idx(p_pe,1)/), count=(/dnparts/))
       IF (rcode /= 0) &
            CALL error_bi('Cannot read variable '//TRIM(uppercase(init_vertcoorname)),substr)
       LEV(dnparts+1:dnparts_max) = mdi 
#endif
#ifdef ECHAM5
       SELECT CASE (TRIM(init_vertcoorname))
          CASE ('press')
             rcode=nf90_inq_varid(ncid,'PRESS_INIT',varid)
          CASE ('theta')
             rcode=nf90_inq_varid(ncid,'THETA',varid)
          CASE ('zeta')
             rcode=nf90_inq_varid(ncid,'ZETA',varid)
          CASE DEFAULT
             call error_bi ("Name of vertical coordinate in INIT/POS file not correct!!!",substr)
          END SELECT
       IF (rcode /= 0) CALL error_bi('Cannot find variable  '//TRIM(init_vertcoorname),substr)
       rcode = nf90_get_var(ncid,varid,LEV(1:dnparts),start=(/idx(p_pe,1)/), count=(/dnparts/))
       IF (rcode /= 0) CALL error_bi('Cannot read variable  '//TRIM(init_vertcoorname),substr)
       LEV(dnparts+1:dnparts_max) = mdi 
#endif
       ! read times
       !write (*,*) 'read time_init'
       rcode = nf90_inq_varid(ncid,'TIME_INIT',varid)
       IF (rcode /= 0) CALL error_bi('Cannot find variable TIME_INIT ',substr)
       rcode = nf90_get_var(ncid,varid,JULSEC(1:dnparts),start=(/idx(p_pe,1)/), count=(/dnparts/))
       IF (rcode /= 0) CALL error_bi('Cannot read variable TIME_INIT ',substr)
       JULSEC(dnparts+1:dnparts_max) = mdi 

       !--------------------------------------------------------------
       ! Store initial positions for MIX(_OLD_MIX) / CHEM(_OLD):
       !--------------------------------------------------------------
       LAT_OLD_MIX   = LAT
       LON_OLD_MIX   = LON
       LAT_OLD       = LAT
       LON_OLD       = LON
       LEV_OLD       = LEV
       JULTIME_OLD   = JULTIME

#ifdef ECHAM5
       CALL clams_locate_parcel_cells(2) ! op_sb_201090801
#endif
       !--------------------------------------------------------------
       ! read chemical species from init file
       !--------------------------------------------------------------

       allocate(field_init(dnparts))
       field_init = mdi

       ! Read species
       do i=1, nspec-ntrac_cl
#ifdef ECHAM5
          IF (init_h2o_emac) THEN
             IF (TRIM(SPECARR(i)%name)=='H2O')     CYCLE
             IF (TRIM(SPECARR(i)%name)=='H2O_100') CYCLE
          END IF
          IF (TRIM(SPECARR(i)%name)=='IWC')     CYCLE
          IF (TRIM(SPECARR(i)%name)=='IWC_100') CYCLE
          IF (TRIM(SPECARR(i)%name)=='CLWC')    CYCLE
#endif
         
          if (p_pe==0) write (*,*) 'read variable in clams_global_start  ', &
               trim(SPECARR(i)%name)
         

!!!!!             
!          SPECARR(i)%values = mdi
          SPECARR(i)%values = 0.
          rcode = nf90_inq_varid(ncid,trim(SPECARR(i)%name),varid)
          if (rcode/=nf90_noerr) then
             if (trim(SPECARR(i)%name)=='H2O') then 
                rcode = nf90_inq_varid(ncid,'SH',varid)
                if (rcode/=nf90_noerr) init_sh = .TRUE.
             elseif (trim(SPECARR(i)%name)=='H2O_100') then 
                rcode = nf90_inq_varid(ncid,'SH_100',varid)
                if (rcode/=nf90_noerr) init_sh = .TRUE.
             endif
          endif
          if (rcode/=nf90_noerr) then
             if (p_pe==0) then
                write (*,*) 'WARNING: variable ',trim(SPECARR(i)%name), &
                     ' not in initfile and set to 0. !!!'
             endif
             CYCLE
          endif
          
          rcode = nf90_get_var(ncid,varid,field_init, &
               start=(/idx(p_pe,1)/), count=(/dnparts/))
          if (rcode/=nf90_noerr) &
               call error_bi(substr,"Cannot read variable "//trim(SPECARR(i)%name))
          
          if (trim(SPECARR(i)%name)/='H2O' .AND. &
               trim(SPECARR(i)%name)/='H2O_100') then 
             SPECARR(i)%values(1:dnparts) = field_init
          else
             if (init_sh) THEN 
                SPECARR(i)%values(1:dnparts) = field_init*28.9644/18.015 
             else
                SPECARR(i)%values(1:dnparts) = field_init
             end if
          end if
        
       enddo

       deallocate(field_init)

!!!!! ???       
          ! ! Read PV_INIT (if available)
          ! PV = mdi
          ! rcode = nf90_inq_varid(ncid,'PV_INIT',varid)0
          ! if (rcode==0) then
          !    rcode = nf90_get_var(ncid,varid,PV(1:dnparts),start=(/idx(p_pe,1)/), count=(/dnparts/))
          ! else 
          !    print *, 'PV_INIT is not available (set on 0.0)'
          !    PV(:)=0.0
          ! endif

       ! Close init file
       rcode=nf90_close(ncid)
       
    end if !lstart

    !*******************************************************
    ! Set switches for submodel events
    !*******************************************************
    IF (USE_CLAMSMIX) THEN
       lmixevent = event_state(mixevent, current_date)
    ELSE
       lmixevent = .FALSE.
    END IF
    
    IF (USE_CLAMSBMIX) THEN
       lbmixevent = event_state(bmixevent, current_date)
    ELSE
       lbmixevent = .FALSE.
    END IF

    IF (USE_CLAMSCIRRUS) THEN
       lcirrusevent = event_state(cirrusevent, current_date)
    ELSE
       lcirrusevent = .FALSE.
    END IF

    lchemevent = .FALSE. ! op_pj_20170110
!!#D clamschem +
    IF (USE_CLAMSCHEM) THEN
       if (timestep_chem == delta_time) then
          lchemevent = .TRUE.
       else
          lchemevent = event_state(chemevent, current_date)
       endif
    ELSE
       lchemevent = .FALSE.
    END IF
!!#D clamschem -

    IF (USE_CLAMSSEDI) THEN
       if (timestep_sedi == delta_time) then
          lsedievent = .TRUE.
       else
          lsedievent = event_state(sedievent, current_date)
       endif
    ELSE
       lsedievent = .FALSE.
    END IF

    IF (USE_CLAMSTRACER) THEN
       if (timestep_tracer*3600 == delta_time) then
          ltracerevent = .TRUE.
       else
          ltracerevent = event_state(tracerevent, current_date)
       endif
    ELSE
       ltracerevent = .FALSE.
    END IF

    IF (USE_CLAMSDEEPCONV) THEN
       if (timestep_deepconv*3600 == delta_time) then
          ldeepconvevent = .TRUE.
       else
          ldeepconvevent = event_state(deepconvevent, current_date)
       endif
    ELSE
       ldeepconvevent = .FALSE.
    END IF


!!!!!
    lclamsoutevent = event_state(clamsoutevent, current_date)

    !*******************************************************
    ! Reset mean age tracer in perpetuum run:
    !*******************************************************
    IF (lperpetuum) THEN
       IF (MONTH==1 .AND. DAY==1 .AND. HOUR==0 .AND. MINUTE==0) THEN
          DO i=1, nspec
             IF (trim(SPECARR(i)%name)/='BA') THEN
                CYCLE
             ELSE
                DO ipart=1,dnparts_max
                   IF (ABS((SPECARR(i)%values(ipart)-mdi)/mdi)>eps) THEN
                      SPECARR(i)%values(ipart) = SPECARR(i)%values(ipart)-0.365
                   END IF
                END DO
             END IF
          END DO
       END IF
    END IF


  END SUBROUTINE clams_global_start

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clams_global_end(flag)

    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_blather_bi,       ONLY: error_bi
#ifdef ECHAM5
    USE messy_main_transform_bi,     ONLY: trp_gpdc_gpgl
    USE messy_main_grid_def_mem_bi,  ONLY: nlev, ngl
    USE messy_main_tracer_mem_bi,    ONLY: xt_c
#endif
 
    ! SMCL: MESSy 
    USE messy_main_channel,            ONLY: set_channel_output, get_channel_object_dimvar
    USE messy_main_timer,              ONLY: lfirst_cycle, lstart, l_rerun, &
                                             YEAR_START, MONTH_START, DAY_START, &
                                             HOUR_START, MINUTE_START, SECOND_START
    USE messy_main_switch,             ONLY: USE_CLAMSMIX
    USE messy_main_tools,              ONLY: PTR_1D_ARRAY
    USE messy_main_constants_mem,      ONLY: STRLEN_ULONG

    ! SMCL: CLaMS
    USE messy_clams_global,    ONLY: DP, mdi, &
                                     lclamsoutevent, nparams, nparams_old, &
                                     PARAM, PARAM_OLD, paramnames, paramnames_old, &
                                     fut_year, fut_month, fut_day, fut_sec, &
                                     irdday, irdsec, dates30, &
                                     lcoupled, init_vertcoorname, ldiagout, &
                                     clams_gridding, &
                                     n_cltr, SPECARR, & ! op_sb_20191010
#ifndef ECHAM5
                                     timetype
#else    
                                     timetype, &
                                     nx, ny, nz, longrid, latgrid, levelgrid, &
                                     longrid_rad, latgrid_rad, asc_lat, &
                                     E5_LAT, E5_LON, E5_LEVEL, &
                                     PREDATA, FUTDATA, &
                                     UDT, VDT, WDT, LEVELDT, &
                                     UFUT, VFUT, WFUT, LEVELFUT
#endif                                     
    USE messy_clamsmix_global,         ONLY: switch_mixing, vert_mix_param
    USE messy_clams_tools_utils,       ONLY: uppercase, bubble_sort, str_pos
    USE messy_clams_tools_interpolreg, ONLY: interpolate_param
    USE messy_clams_tools_dateconv,    ONLY: incdat_interface

    IMPLICIT NONE

   
    INTEGER, INTENT(IN) :: flag   ! op_sb_20200402
    
    CHARACTER(LEN=*), PARAMETER :: substr = 'clams_global_end'
    REAL(DP)       :: pi
    INTEGER        :: status
    INTEGER        :: i, j, k, iparam, jt, i1
    TYPE(timetype) :: itime, ipasttime
    
#ifdef ECHAM5    
    TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER :: dva  => NULL()
    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER :: units => NULL()
#endif    

    pi=4.*atan(1.) 
   SELECT CASE(flag)
    CASE(1)
#ifdef ECHAM5 
! op_sb_20191105+
    ALLOCATE(KHPBL(nlon,ngl))   
    KHPBL(:,:)=nlev-2                                                    
! op_sb_20191105-
    !--------------------------------------------------------------
    ! ECHAM WIND
    !
    ! => UDT, VDT, WDT, LEVELDT
    ! => UFUT, VFUT, WFUT, LEVELFUT
    ! => PREDATA, FUTDATA
    !--------------------------------------------------------------
    
    ! Collect global fields:
    CALL trp_gpdc_gpgl(1, E5_LAT3D_D , E5_LAT3D_G)
    CALL trp_gpdc_gpgl(1, E5_LON3D_D , E5_LON3D_G)
    CALL trp_gpdc_gpgl(1, E5_UWIND_D , E5_UWIND_G)
    CALL trp_gpdc_gpgl(1, E5_VWIND_D , E5_VWIND_G)
    CALL trp_gpdc_gpgl(1, E5_TEMP_D  , E5_TEMP_G )
    CALL trp_gpdc_gpgl(1, E5_PRESS_D , E5_PRESS_G)
    CALL trp_gpdc_gpgl(1, E5_THETA_D, E5_THETA_G)
    CALL trp_gpdc_gpgl(1, E5_PSURF_D, E5_PSURF_G)
    CALL trp_gpdc_gpgl(1, E5_OMEGA_D, E5_OMEGA_G)
    CALL trp_gpdc_gpgl(1, E5_TTE_D, E5_TTE_G)
    CALL trp_gpdc_gpgl(1, E5_ALPSTE_D, E5_ALPSTE_G)
    CALL trp_gpdc_gpgl(1, E5_ZETA_D, E5_ZETA_G)
    CALL trp_gpdc_gpgl(1, E5_ZETADOT_D, E5_ZETADOT_G)
    CALL trp_gpdc_gpgl(1, E5_SIGMA_D, E5_SIGMA_G)
    CALL trp_gpdc_gpgl(1, E5_SIGMADOT_D, E5_SIGMADOT_G)
    CALL trp_gpdc_gpgl(1, E5_THETADOT_D, E5_THETADOT_G)
    CALL trp_gpdc_gpgl(1, E5_THETADOT2_D, E5_THETADOT2_G)
    CALL trp_gpdc_gpgl(1, E5_HEATING_D, E5_HEATING_G)
    CALL trp_gpdc_gpgl(1, E5_PSDOT_D, E5_PSDOT_G)
    CALL trp_gpdc_gpgl(1, E5_SH_D, E5_SH_G)
    CALL trp_gpdc_gpgl(1, E5_IWC_D, E5_IWC_G)
    CALL trp_gpdc_gpgl(1, E5_CLWC_D, E5_CLWC_G)
    ! op_sb_20191105+
     CALL trp_gpdc_gpgl(1, pblh_i, khpbl)    
    ! decomposition -> global field (grid mass)
     CALL trp_gpdc_gpgl(1, grmass, pmbox)
    ! op_sb_20191105-
   
! op_pj_20180614+
! Note: The transposition from decomposed to global is time (MPI all-to-all
!       communication) and memory intensive (global target fields on all 
!       tasks). Therefore it is more efficient to perform local calculations
!       with the decomposed fields where possible (e.g., the summation below)
!       and only transpose the required result.
!       Yet, with the new approach to extract the diabatic heating rate
!       this all becomes obsolete ...
!
!!$    CALL trp_gpdc_gpgl(1, E5_HEAT_VDIFF_D, E5_HEAT_VDIFF_G)
!!$    CALL trp_gpdc_gpgl(1, E5_HEAT_CLOUDCONV_D, E5_HEAT_CLOUDCONV_G)
!!$    CALL trp_gpdc_gpgl(1, E5_HEAT_RHEAT_D, E5_HEAT_RHEAT_G)
!!$    CALL trp_gpdc_gpgl(1, E5_HEAT_SSO_D, E5_HEAT_SSO_G)
! op_pj_20180614-

    nx = SIZE(E5_UWIND_G,1)
    ny = SIZE(E5_UWIND_G,3)
  
    CALL get_channel_object_dimvar(status, 'g2a', 'um1', dva, units)

    E5_LON   => dva(1)%ptr(:)
    E5_LAT   => dva(3)%ptr(:)
    E5_LEVEL => dva(2)%ptr(:)

    ! 3D LAT&LON ARRAYS for ECHAM    
    DO k=1,nlev
       DO i=1,nx
          E5_LAT3D_G(i,k,:) = E5_LAT(:)
       END DO
       DO j=1,ny
          E5_LON3D_G(:,k,j) = E5_LON(:)
       END DO
    END DO

!!$    IF (lfirst_cycle) THEN  !Allocate target-arrays: 
!!$       ALLOCATE (UDT_T (nx, ny, nlev))
!!$       ALLOCATE (VDT_T (nx, ny, nlev))
!!$       ALLOCATE (WDT_T (nx, ny, nlev))
!!$    END IF

! op_pj_20180614+
!!$    ! Calculate total heating rate Q:
!!$    E5_HEATING_G = E5_HEAT_VDIFF_G + E5_HEAT_CLOUDCONV_G + &
!!$                   E5_HEAT_RHEAT_G + E5_HEAT_SSO_G
! op_pj_20180614-

    ! Calculate theta_dot:
    CALL calculate_thetadot_echam
    ! Calculate zeta & zetadot:
    CALL calculate_zeta_zetadot_echam(nx, ny, nlev)
    WHERE (E5_ZETA_G .LT. 0.)
       E5_ZETA_G = 0.
    END WHERE

    ! Switch dimensions for TRAJ input fields:
    ! Change units of pressure (Pa -> hPa) and pressure tendency (Pa/s -> hPa/s)

    DO k=1,nlev                                  ! Loop over levels
       UDT(:,:,nlev-k+1) = E5_UWIND_G(:,k,:) / COS(E5_LAT3D_G(:,k,:)*pi/180.)
       VDT(:,:,nlev-k+1) = E5_VWIND_G(:,k,:) / COS(E5_LAT3D_G(:,k,:)*pi/180.)
    END DO

    SELECT CASE (TRIM(init_vertcoorname))
    CASE ('press')
    DO k=1,nlev
       WDT(:,:,nlev-k+1) = E5_OMEGA_G(:,k,:) / 100. 
       leveldt(:,:,nlev-k+1) = E5_PRESS_G(:,k,:)/ 100.
       levelfut(:,:,nlev-k+1) = E5_PRESS_G(:,k,:)/ 100.
    END DO
    CASE ('theta')
    DO k=1,nlev
       WDT(:,:,nlev-k+1) = E5_THETADOT_G(:,k,:)/86400.
       leveldt(:,:,nlev-k+1) = E5_THETA_G(:,k,:)
       levelfut(:,:,nlev-k+1) = E5_THETA_G(:,k,:)
    END DO
    CASE ('zeta')
    DO k=1,nlev
       WDT(:,:,nlev-k+1) = E5_ZETADOT_G(:,k,:)/86400.
       leveldt(:,:,nlev-k+1) = E5_ZETA_G(:,k,:)
       levelfut(:,:,nlev-k+1) = E5_ZETA_G(:,k,:)
    END DO
    CASE DEFAULT
       call error_bi ("Name of vertical coordinate in INIT/POS file not correct!!!",substr)
    END SELECT

    DO iparam = 1, nparams
       SELECT CASE (PREDATA(iparam)%name)
       CASE ('TEMP')
          DO k=1,nlev
          PREDATA(iparam)%values(:,:,nlev-k+1) = E5_TEMP_G(:,k,:) 
          FUTDATA(iparam)%values(:,:,nlev-k+1) = E5_TEMP_G(:,k,:) 
          END DO
       CASE ('PRESS')
          DO k=1,nlev
          PREDATA(iparam)%values(:,:,nlev-k+1) = E5_PRESS_G(:,k,:)/ 100.
          FUTDATA(iparam)%values(:,:,nlev-k+1) = E5_PRESS_G(:,k,:)/ 100.
          END DO
       CASE ('H2O')
          DO k=1,nlev
          PREDATA(iparam)%values(:,:,nlev-k+1) = E5_SH_G(:,k,:) &
               *28.9644/18.015
          FUTDATA(iparam)%values(:,:,nlev-k+1) = E5_SH_G(:,k,:) &
               *28.9644/18.015
          END DO
       CASE ('IWC')
          DO k=1,nlev
          PREDATA(iparam)%values(:,:,nlev-k+1) = E5_IWC_G(:,k,:)
          FUTDATA(iparam)%values(:,:,nlev-k+1) = E5_IWC_G(:,k,:)
          END DO
       CASE ('CLWC')
          DO k=1,nlev
          PREDATA(iparam)%values(:,:,nlev-k+1) = E5_CLWC_G(:,k,:)
          FUTDATA(iparam)%values(:,:,nlev-k+1) = E5_CLWC_G(:,k,:)
          END DO
       CASE ('ZETA')
          DO k=1,nlev
          PREDATA(iparam)%values(:,:,nlev-k+1) = E5_ZETA_G(:,k,:)
          FUTDATA(iparam)%values(:,:,nlev-k+1) = E5_ZETA_G(:,k,:)
          END DO
       CASE ('ZETA_DOT')
          DO k=1,nlev
          PREDATA(iparam)%values(:,:,nlev-k+1) = E5_ZETADOT_G(:,k,:)
          FUTDATA(iparam)%values(:,:,nlev-k+1) = E5_ZETADOT_G(:,k,:)
          END DO
       CASE ('U')
          DO k=1,nlev
          PREDATA(iparam)%values(:,:,nlev-k+1) = E5_UWIND_G(:,k,:)
          FUTDATA(iparam)%values(:,:,nlev-k+1) = E5_UWIND_G(:,k,:)
          END DO
       CASE ('V')
          DO k=1,nlev
          PREDATA(iparam)%values(:,:,nlev-k+1) = E5_VWIND_G(:,k,:)
          FUTDATA(iparam)%values(:,:,nlev-k+1) = E5_VWIND_G(:,k,:)
          END DO
       CASE ('THETA')
          DO k=1,nlev
          PREDATA(iparam)%values(:,:,nlev-k+1) = E5_THETA_G(:,k,:)
          FUTDATA(iparam)%values(:,:,nlev-k+1) = E5_THETA_G(:,k,:)
          END DO
       CASE ('PV')
          PREDATA(iparam)%values(:,:,:) = 1.
          FUTDATA(iparam)%values(:,:,:) = 1.
       CASE DEFAULT
          call error_bi (TRIM(PREDATA(iparam)%name)//' not initialized !!!',substr)
       END SELECT
    END DO

    ! PREDATA=FUTDATA and UDT=UFUT etc. for EMAC INPUT 
    ! (no interpolation in time in TRAJ) 
!!$    UDT  => UDT_T
!!$    VDT  => VDT_T
!!$    WDT  => WDT_T
!!$    UFUT => UDT_T
!!$    VFUT => VDT_T
!!$    WFUT => WDT_T
    UFUT = UDT
    VFUT = VDT
    WFUT = WDT

    CALL trp_gpdc_gpgl(-1, E5_LAT3D_D, E5_LAT3D_G)
    CALL trp_gpdc_gpgl(-1, E5_LON3D_D, E5_LON3D_G)
    CALL trp_gpdc_gpgl(-1, E5_ZETA_D, E5_ZETA_G)
    CALL trp_gpdc_gpgl(-1, E5_ZETADOT_D, E5_ZETADOT_G)
    CALL trp_gpdc_gpgl(-1, E5_SIGMA_D, E5_SIGMA_G)
    CALL trp_gpdc_gpgl(-1, E5_SIGMADOT_D, E5_SIGMADOT_G)
    CALL trp_gpdc_gpgl(-1, E5_THETADOT_D, E5_THETADOT_G)
    CALL trp_gpdc_gpgl(-1, E5_THETADOT2_D, E5_THETADOT2_G)
! op_pj_20180614+: Note: This back-transposition from global to decomposed
!                        is most probably not needes, since E5_HEATING_G
!                        was not changed ...?
    CALL trp_gpdc_gpgl(-1, E5_HEATING_D, E5_HEATING_G)
! op_pj_20180614-
    CALL trp_gpdc_gpgl(-1, E5_PSDOT_D, E5_PSDOT_G)
#endif


#ifdef ECHAM5
    
    !--------------------------------------------------------------
    ! Set levelgrid, latgrid & longrid
    ! (MBM-CLaMS: in clams_setup)
    !--------------------------------------------------------------
    IF (lfirst_cycle) THEN
       nx = SIZE(E5_LON)
       ny = SIZE(E5_LAT)
       nz = SIZE(E5_LEVEL)

       ALLOCATE(levelgrid(nz)) 
       DO i=1,nz
          levelgrid(i) = E5_LEVEL(i)
       END DO

       allocate (longrid(nx))
       allocate (latgrid(ny))
       longrid(1:nx) = E5_LON
       latgrid(1:ny) = E5_LAT

!!!!!  bubble_sort necessary ???
       CALL bubble_sort (levelgrid,nz)

       if (p_pe==0) then
          WRITE (*,*)
          WRITE (*,'(A,I4)') 'Number of longitudes: ', nx
          WRITE (*,'(A,I4)') 'Number of latitudes: ', ny
          WRITE (*,'(A,I4)') 'Number of levels: ', nz
          WRITE (*,*)
          if (ldiagout) then
             write (*,*) 'longrid=',longrid
             write (*,*) 'latgrid=',latgrid
             write (*,*) 'levelgrid=',levelgrid
             WRITE (*,*)
          endif
       endif

       allocate (longrid_rad(nx+1))
       allocate (latgrid_rad(0:ny+1))

       longrid_rad(1:nx) = longrid(1:nx)/180. * pi
       latgrid_rad(1:ny) = latgrid(1:ny)/180. * pi

       longrid_rad(nx+1) = (longrid(1)+360.)/180. * pi
       if (latgrid(2)>latgrid(1)) then
          asc_lat = .true.
          latgrid_rad(0)    = -0.5 * pi
          latgrid_rad(ny+1) = 0.5 * pi
       else
          asc_lat = .false.
          latgrid_rad(0)    = 0.5 * pi
          latgrid_rad(ny+1) = -0.5 * pi
       endif
    END IF
#endif

    !--------------------------------------------------------------
    ! Interpolate parameters for start time
    !--------------------------------------------------------------
    if (lstart) then

       IF (p_pe==0) write(*,*) uppercase(substr), ': in lstart'


       ! current time
       itime%year    = YEAR_START
       itime%month   = MONTH_START
       itime%day     = DAY_START
       itime%sec     = HOUR_START*3600+MINUTE_START*60+SECOND_START

       ! Calculate time of previous windfile
       ipasttime%sec = fut_sec
       ipasttime%day = fut_day
       ipasttime%month = fut_month
       ipasttime%year = fut_year

       CALL incdat_interface( ipasttime%sec,ipasttime%day,ipasttime%month,ipasttime%year, &
            -irdsec,-irdday,0,0, dates30)     

       DO i = 1, nparams

          PARAM(i)%values     = mdi

          CALL interpolate_param (LON, LAT, LEV, itime, ipasttime, &
               i, paramnames(i), PARAM(i)%values, lcoupled)

       ENDDO

       DO i = 1, nparams_old
          PARAM_OLD(i)%values = PARAM(str_pos(nparams,paramnames,paramnames_old(i)))%values
          !write (*,*) 'iparam_old, iparam:', i, str_pos(nparams,paramnames,paramnames_old(i))
       ENDDO

       ! if MIX is used and vertical mixing is switched on:
       !    initialize THETA_OLD_MIX and BVF_OLD_MIX
       ! THETA and BVF must be specified as parameter (saved on PARAM)
       IF (USE_CLAMSMIX .and. switch_mixing==2) THEN 
          THETA_OLD_MIX  = PARAM(str_pos(nparams,paramnames,'THETA'))%values
          if (uppercase(vert_mix_param)=='WET') THEN
             BVF_OLD_MIX = PARAM(str_pos(nparams,paramnames,'BVF_WET'))%values
          else
             BVF_OLD_MIX = PARAM(str_pos(nparams,paramnames,'BVF'))%values
          endif
       ENDIF

    endif

!!$    if (lclamsoutevent) then
!!$       CALL set_channel_output(status, 'clams', .TRUE.)
!!$       CALL channel_halt(substr, status)
!!$       
!!$       ! wind data output: 
!!$       CALL set_channel_output(status, 'winddata', .TRUE.)
!!$       CALL channel_halt(substr, status)
!!$       
!!$!!!!! => in  bml/clams_main.f90      
!!$!       CALL messy_write_output
!!$    endif

!!!!! => in  bml/clams_main.f90      
!    IF (l_rerun) CALL messy_write_restart

#ifndef NOMPI
#ifdef ECHAM5
    IF (clams_gridding) CALL clams_tracer_update ! ju_ec_20180628
!!$    CALL clams_update_celldist ! op_sb_2019071
! op_sb_20191010+
  ! Update MESSy-TRACER Lagrangian representation
!     DO i1 = 1,n_cltr
!      jt = CLTR_IDX(i1)
!      xt_c(:,1,jt,1) = SPECARR(CLTR_SPECARR_IND(i1))%values
!     ENDDO
! op_sb_20191010-
    CALL clams_locate_parcel_cells(2) ! op_sb_20190801
#endif
#endif
   ! op_sb_20200402+
   CASE(2)
#ifdef ECHAM5
    CALL clams_locate_parcel_cells(2)
#endif
   END SELECT
   ! op_sb_20200402-
  END SUBROUTINE clams_global_end
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clams_free_memory

    USE messy_clams_global, ONLY: SPECARR, PARAM, PARAM_OLD, &
                                  PREDATA, FUTDATA, &
                                  corr_thetadot, corr_3d, &
                                  time_corr, theta_corr, lat_corr, &
                                  thetadot_corr, thetadot_corr3d, &
                                  levelgrid, latgrid, longrid, &
                                  latgrid_rad, longrid_rad
    USE messy_main_mpi_bi,  ONLY: p_pe

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clams_free_memory'

    IF (p_pe==0) write(*,*) substr

    if (associated(SPECARR)) deallocate (SPECARR)

    if (associated(PREDATA)) DEALLOCATE (PREDATA)
    if (associated(FUTDATA)) DEALLOCATE (FUTDATA)

    if (associated(PARAM))     DEALLOCATE (PARAM)
    if (associated(PARAM_OLD)) DEALLOCATE (PARAM_OLD)

    if (associated(levelgrid))   DEALLOCATE (levelgrid)
    if (associated(latgrid))     DEALLOCATE (latgrid)
    if (associated(longrid))     DEALLOCATE (longrid)
    if (associated(latgrid_rad)) DEALLOCATE (latgrid_rad)
    if (associated(longrid_rad)) DEALLOCATE (longrid_rad)
   
#ifdef ECHAM
    if (associated(pc_d)) then
       DEALLOCATE (pc_d)
       NULLIFY(pc_d)
    endif
    if (associated(pc_g)) then
       DEALLOCATE (pc_g)
       NULLIFY(pc_g)
    endif
     if (associated(clgp_data)) then
       DEALLOCATE (clgp_data)
       NULLIFY(clgp_data)
    endif  
    if (associated(pressi_init)) then
       DEALLOCATE (pressi_init)
       NULLIFY(pressi_ini)
    endif 
#endif
    ! op_sb_20191105+
     if (associated(khpbl)) then
        DEALLOCATE(khpbl) 
        NULLIFY(khpbl)
     endif 
! op_sb_20191105-
#ifdef MBM_CLAMS    
    if (corr_thetadot) then
       deallocate (time_corr)
       deallocate (theta_corr)
       if (corr_3d) then
          deallocate (lat_corr)
          deallocate (thetadot_corr3d)
       else
          deallocate (thetadot_corr)
       endif
    endif
#endif  

!!$#ifdef ECHAM5
!!$    DEALLOCATE (UDT_T)
!!$    DEALLOCATE (VDT_T)
!!$    DEALLOCATE (WDT_T)
!!$#endif  

  END SUBROUTINE clams_free_memory


!--------------------------------------------------------------------

!--------------------------------------------------------------------
!
! ju_ec_20190212 
! COMMENT REGARDING COUPLING BETWEEN CLAMS AND EMAC
! Edward Charlesworth
! 2019.02.12
!
! These comments are written to help you, the reader, understand the
! coupling between EMAC and CLAMS. I struggled for a long time during
! my PhD to identify errors in the EMAC->CLaMS coupling, and would
! have been helped greatly by good commenting in the related
! routines. Therefore I am now trying my best to comment this work
! for future coders. So, before you read these comments, I want to 
! ask you to, please, write extensive comments in your own code. That
! will save time and energy for your collegaues, for your successors,
! and for you (because you'll have a reminder of what you've done!).
!
! What should you commend? Comment everything! Write comments on the 
! variables and their units, write a table of them at the start of 
! routines, give your variables obvious and useful names, describe 
! what your routines do,  explain the changes you think could or 
! should be made, and write whatever else you think about your code. 
! It's right thing to do, and a healthy, productive modeling group 
! will be greatly helped by that work.
!
! So thank you for commenting, and good luck modeling!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! What does the coupling involve? There are two components to the 
! coupling:
! 
!   a. The transport coupling was the primary work of Charlotte Hoppe 
!      during her PhD in the early half of the 2010's. Nicole Thomas
!      and Patrick Joeckel also worked closely on this. The transport
!      coupling takes EMAC variables for wind and heating rates and
!      applies them as the CLaMS advection variables.
!
!   b. The radiative-feedback coupling was the primary work of my PhD
!      while Nicole Thomas and Patrick Joeckel were also closely
!      involved. The radiative-feedback coupling applies the chemical 
!      fields from CLaMS to EMAC radiation.
!
! These components are described in the following comment sections.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! WHAT THE TRANSPORT COUPLING IS
!
! The transport coupling involves acquisition of the relevant EMAC
! variables via get_channel_object calls, transformation of the 
! variables via trp_gpdc_gpgl (in the routines above), and then
! computation of the relavant diabatic transport variables (zeta and
! its tendency). I think these variables are also transformed into
! the decomposed versions via trp_gpdc_gpgl? Eventually, these 
! variables are used for the trajectory routines, but I don't know
! precisely how that happens.
!
! The transport coupling subroutines (that I know of) are:
!
!   1. calculate_thetadot_echam: This calculates the theta tendency,
!      which is later used in the zeta tendency calculation.
!
!   2. calculate_zeta_zetadot_echam: This routine calculates the
!      necessary zeta and zeta tendency variables.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! WHAT THE RADIATIVE-FEEDBACK COUPLING IS
!
! The radiative-feedback coupling involves two major components -
! Lagrangian-to-grid transformation and grid-to-Lagrangian 
! transformation - and a major new MESSy infrastructure integration
! into MESSy-CLaMS. The infrastructure integrated is MESSy-TRACER,
! for which two tracer sets were created. Those are "cl" and "clgp".
! Both of these have channels with "tracer_" preceding their names.
! The first tracer set contains the Lagrangian representation of
! CLaMS (it's redundant with the normal CLaMS channel). The second
! contains the gridded data.
!
! The clgp tracer set contains CLaMS data gridded to the EMAC grid. 
! The process of this gridding is "binning", as opposed to the normal
! CLaMS approach of interpolation. Binning is faster and proceeds as:
! the EMAC grid cell which contains each parcel is identified, the 
! clgp 3d field (for each tracer) is initialized with zeros, and a DO
! loop runs over the parcels, adding their data to their respective
! grid cells (divided by the number of parcels in the cell). A
! disadvantage of the binning process is that some grid cells in
! every step will be empty of parcels, and so their data is filled
! either with negative values or with the values from the EMAC "fill 
! channel object". For the coupling, a fill object is required to 
! prevent  these empty cells. The fill object will also be used for 
! one more process. Oh yeah, this fill object is an object in some
! channel. It could be any channel and any object within that channel
! as long as it's a 3d field.
!
! Precisely, the coupling cannot reasonably occur in the entire
! model domain, so a "clams region" must be defined, and outside of
! this region the clams parcels must be filled with the fill object
! data. That's right, the parcels themselves are filled! Effectively,
! this non-clams region is the boundary condition for clams, then.
! The reason for this is to prevent sharp gradients in the clgp field
! but this also has the advantage of providing reasonable 
! tropospheric values of some tracers to the clgp field, as long as
! the fill object is reasonable in the troposphere.
!
! Right, so what are the subroutines involved?
!
!   1. clams_new_tracer: Adds tracers to the grid and parcel tracer 
!      sets. Maybe more coding could be done there in the future?
!
!!$!   2. clams_init_tracer: Just initializes the tracer sets with 
!!$!      zeros.
!
!   3. clams_update_tracer: This handles the boring parts of updating
!      the clgp fields and the MESSy-TRACER clams parcel field. Calls
!      to more meaty routines, parrallelization transformations, and
!      so on go here.
!
!   4. clams_locate_parcel_cells: Assigns location indexes to the
!      parcels. These are stored in the variable POS.
!
!   5. clams_lg2gp: Bins the clams parcel data onto the clgp grid. 
!
!   6. clams_gp2lg: Replaces parcel data outside for parcels outside
!      the clams region with the fill channel object data for the
!      relevant tracer.
!
!
!--------------------------------------------------------------------


!--------------------------------------------------------------------
  SUBROUTINE calculate_thetadot_echam
  ! Calculate vertical velocity theta_dot from ECHAM input fields

  !
  ! ju_ec_20190214
  !
  !    |    Variable    |   units   |    computed in / obtained from      |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_HEATING_G    |   K / s   | from object tte of scnbuf           |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_TEMP_G       |     K     | object 'tm1' of channel 'ECHAM5'    |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_THETA_G      |     K     | object 'tpot' of channel 'ECHAM5'   |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_THETADOT_G   |   K / s   | here                                |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_THETADOT2_G  |   K / s   | here, I don't know what this is for |
  !    |                |           | and I'm pretty sure it's pointless  |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |                |           |                                     |
  !

    IMPLICIT NONE

    ! Heatingrates: 
    E5_THETADOT_G = E5_HEATING_G * E5_THETA_G / E5_TEMP_G

    ! Pressure tendency and temperature tendency:
    E5_THETADOT2_G = E5_TTE_G * (100000./E5_PRESS_G)**0.286 - &
                    E5_THETA_G * E5_OMEGA_G / (0.286 * E5_PRESS_G)
    ! ju_ec_20190214: I haven't documented the variables for E5_THEATDOT2
    ! because I don't think it's used for anything. I guess it was an
    ! alternative computation that was later set aside?

  END SUBROUTINE calculate_thetadot_echam
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE calculate_zeta_zetadot_echam(nx, ny, nlev)
  ! Calculate vertical coordinate zeta from ECHAM input fields

  !
  ! ju_ec_20190129
  !
  ! If you modify this code, please update this comment. Thereby, you will
  ! greatly assist all persons working on or interested in this who follow you.
  ! Thank you very much.
  !
  !    |    Variable    |   units   |    computed in / obtained from      |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_ALPSTE_G     |   1 / s   | object 'alpste' of channel 'scnbuf' |
  !    |(tendency of log|           |                                     |
  !    | of surface P)  |           |                                     |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_OMEGA_G      |  Pa / s   | object 'vervel' of channel 'scnbuf' |
  !    |(vert. velocity |           |                                     |
  !    |  in pressure)  |           |                                     |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_PSDOT_5      |  Pa / s   | this routine                        |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_PSURF_G      |    Pa     | object 'aps' of channel 'g3b'       |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_PRESS_G      |    Pa     | object 'press' of channel 'ECHAM5'  |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_SIGMA_G      | unitless  | this routine                        |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_SIGMADOT_G   |   1 / s   | this routine                        |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_THETA_G      |     K     | object 'tpot' of channel 'ECHAM5'   |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_THETADOT_G   |   K / s   | calculate_thetadot_echam, above     |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_ZETA_G       |     K     | this routine                        |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |E5_ZETADOT_G    |   K / day | this routine                        |
  !    |(vert. velocity |           |                                     |
  !    | used in eclams)|           |                                     |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |dfunc           | unitless  | this routine                        |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |func            | unitless  | this routine                        |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |p_ref           |    hPa    | messy_clams_global, default 300,    |
  !    |                |           | can be defined in clams.nml         |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |delta_time      |  seconds  | a main messy variable, set by the   |
  !    |(global t-step) |           | the namelist timer.nml              |
  ! ---|----------------|-----------|-------------------------------------|---
  !    |                |           |                                     |
  !

    USE messy_clamstraj_global, ONLY: p_ref ! reference pressure level in hPa

    USE messy_main_timer,          ONLY: delta_time

    IMPLICIT NONE
  
    INTEGER, INTENT(IN) :: nx, ny, nlev
    INTEGER             :: i, j, k
    REAL(DP)            :: pi
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: func, dfunc, sigma_ref 

    ! Calculate sigma (= p/p_surface)
    DO k=1,nlev
       E5_SIGMA_G(:,k,:) = E5_PRESS_G(:,k,:) / E5_PSURF_G
    END DO
    
    ! Calculate surface pressure derivative
    E5_PSDOT_G = E5_PSURF_G * ( EXP(E5_ALPSTE_G*delta_time) - 1. )/delta_time

    ! Calculate sigma_dot
    DO k=1,nlev
       E5_SIGMADOT_G(:,k,:) =   E5_OMEGA_G(:,k,:) / E5_PSURF_G &
                            & - E5_PRESS_G(:,k,:) * E5_PSDOT_G / (E5_PSURF_G**2)
    END DO

    ! Calculate function f
    ALLOCATE (func(nx, nlev, ny))
    ALLOCATE (sigma_ref(nx, nlev, ny))
    pi = 4.*ATAN(1.)  

    DO i=1,nx
       DO k=1,nlev
          DO j=1,ny
             sigma_ref(i,k,j)= p_ref*100. / E5_PSURF_G(i,j)
             IF (E5_SIGMA_G(i,k,j) .GT. sigma_ref(i,k,j)) THEN
                func(i,k,j) = SIN( 0.5*pi                       &
                                 & * (1. - E5_SIGMA_G(i,k,j))   &
                                 & / (1. - sigma_ref(i,k,j))    &
                                 & )
             ELSE
                func(i,k,j) = 1.
             END IF
          END DO
       END DO
    END DO

    ! Calculate derivative (d f / d sigma)
    ALLOCATE (dfunc(nx, nlev, ny))
    pi = 4.*ATAN(1.)  
    DO i=1,nx
       DO k=1,nlev
          DO j=1,ny
             IF (E5_SIGMA_G(i,k,j) .GT. sigma_ref(i,k,j)) THEN
                dfunc(i,k,j) = -0.5*pi / (1-sigma_ref(i,k,j))   &
                               & * COS(0.5*pi                   &
                               & * (1. - E5_SIGMA_G(i,k,j))     &
                               & / (1. - sigma_ref(i,k,j))      &
                               & )
             ELSE
                dfunc(i,k,j) = 0.
             END IF
          END DO
       END DO
    END DO
  
    ! Calculate zeta
    DO i=1,nx
       DO k=1,nlev
          DO j=1,ny
             E5_ZETA_G(i,k,j) = E5_THETA_G(i,k,j) * func(i,k,j) 
          END DO
       END DO
    END DO

    ! Calculate zeta_dot
    E5_ZETADOT_G = ( func * E5_THETADOT_G + & 
                   E5_THETA_G * dfunc * E5_SIGMADOT_G) * 86400.

    DEALLOCATE(func)
    DEALLOCATE(dfunc)
    DEALLOCATE(sigma_ref)

  END SUBROUTINE calculate_zeta_zetadot_echam


  ! THESE ARE SOME OLDER CODE SECTIONS FOR THE ABOVE SUBROUTINE
  !DO i=1,nx
  !   DO k=1,nlev
  !      DO j=1,ny
  !         E5_SIGMA_G(i,k,j) = E5_PRESS_G(i,k,j) / E5_PSURF_G(i,j)
  !      END DO
  !   END DO
  !END DO
  !DO i=1,nx
  !   DO k=1,nlev
  !      DO j=1,ny
  !         E5_SIGMADOT_G(i,k,j) = E5_OMEGA_G(i,k,j) / E5_PSURF_G(i,j) -&
  !              & E5_PRESS_G(i,k,j) *E5_PSDOT_G(i,j) / (E5_PSURF_G(i,j)**2)
  !      END DO
  !   END DO
  !END DO
  
#ifdef ECHAM5
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  ! op_pj_20190920: this is not requird at all, beause MESSy tracers
  !                 are initialized with zero anyway
!!$
!!$  !ju_ec_20180627+
!!$  SUBROUTINE clams_init_tracer
!!$
!!$      !
!!$      ! ju_ec_20190214
!!$      ! 
!!$      ! This routine initializes the clams tracer sets with zeros. That's it!
!!$      !
!!$    
!!$  USE messy_main_tracer_mem_bi, ONLY: xt_c, xt_clgp,&
!!$                                      CLTRSTR, CLGPTRSTR !, xtte_c, xtm1_c
!!$  USE messy_main_tracer_bi,     ONLY: tracer_halt
!!$  USE messy_main_tracer,        ONLY: tracer_iniflag
!!$  USE messy_clams_global,       ONLY: n_cltr
!!$
!!$  IMPLICIT NONE
!!$
!!$  ! LOCAL
!!$  INTEGER :: status
!!$  INTEGER :: i, jt
!!$  CHARACTER(LEN=*), PARAMETER :: substr = 'clams_init_tracer'
!!$
!!$  DO i = 1, n_cltr
!!$
!!$    jt = CLTR_IDX(i)
!!$
!!$    xt_c(:,1,jt,1) = 0.0_dp
!!$
!!$    CALL tracer_iniflag(status, CLTRSTR,   jt, .TRUE. )
!!$    CALL tracer_halt( substr, status )
!!$
!!$    xt_clgp(:,:,jt,:) = 0.0_dp
!!$    
!!$    CALL tracer_iniflag(status, CLGPTRSTR, jt, .TRUE. )
!!$    CALL tracer_halt( substr, status )
!!$  
!!$  END DO
!!$  
!!$  END SUBROUTINE clams_init_tracer
!!$!ju_ec_20180627-
!!$!--------------------------------------------------------------------

  !--------------------------------------------------------------------
  SUBROUTINE clams_init_tracer
  
   ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_pe, p_sum
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: main_channel_write_output
    USE messy_main_timer_bi,         ONLY: event_state

    ! SMIL: CLaMS
    USE messy_clamsmix_si,     ONLY: mixevent
    USE messy_clamsbmix_si,    ONLY: bmixevent
    USE messy_clamscirrus_si,  ONLY: cirrusevent

    USE messy_clamssedi_si,    ONLY: sedievent
    USE messy_clamstracer_si,  ONLY: tracerevent

    ! SMCL: MESSy
    USE messy_main_timer,      ONLY: lstart, lresume, current_date, &
                                     YEAR, MONTH, DAY, HOUR, MINUTE,  &
                                     YEAR_START, MONTH_START, DAY_START, &
                                     HOUR_START, MINUTE_START, SECOND_START, &
                                     delta_time, lfirst_cycle
    USE messy_main_switch,     ONLY: USE_CLAMSTRAJ, USE_CLAMSMIX, &
                                     USE_CLAMSBMIX,   &
                                     USE_CLAMSCIRRUS, &
                                     USE_DISSOC, USE_CLAMS, &
                                     USE_CLAMSSEDI, USE_CLAMSTRACER

    USE messy_main_channel,    ONLY: set_channel_output

    ! SMCL: CLaMS
    USE messy_clams_global,    ONLY: DP, prec, mdi, eps, &
                                     initfile, init_vertcoorname, &
                                     nspec, SPECARR,  &
                                     dnparts, dnparts_max, dnparts_co, dnparts_max_co, &
                                     nparts, idx, &
                                     lmixevent, lbmixevent, lcirrusevent, &
                                     lchemevent, lsedievent, ltracerevent, &
                                     lclamsoutevent, &
                                     init_h2o_emac,&
                                     lperpetuum, ldiagout, &
                                     pre_year, pre_month, pre_day, pre_sec, &
                                     pre_year_co, pre_month_co, pre_day_co, pre_sec_co, &
                                     irdday, irdsec, &
                                     fut_year, fut_month, fut_day, fut_sec, &
                                     buffersize

    USE messy_clamssedi_global,    ONLY: timestep_sedi
    USE messy_clamstracer_global,  ONLY: timestep_tracer
    USE messy_clams_read_metdata,  ONLY: clams_read_ecmwf_files
    USE messy_clams_tools_utils,   ONLY: uppercase

    ! op_sb_20200303+
    USE messy_main_grid_def_mem_bi, ONLY: vct, apzero, nvclev,        &
                                          ngl, nlon, nlev
    ! op_sb_20200303-

    USE netcdf

    IMPLICIT NONE

    REAL(PREC), ALLOCATABLE, DIMENSION(:) :: field_init
    REAL(DP),DIMENSION(1)                 :: help_array
    INTEGER :: ncid, rcode, varid
    INTEGER :: i, ipart,j,k
    INTEGER :: status
    INTEGER :: JG,JL,JK

    CHARACTER(LEN=*), PARAMETER :: substr = 'clams_init_tracer'

#ifdef ECHAM5
       IF (.NOT. lresume) then
       IF (p_pe==0) write(*,*) substr
 
       ! Open init file
       rcode = nf90_open(initfile, nf90_nowrite, ncid, buffersize)
       if (rcode /= nf90_noerr) call error_bi(substr,"Cannot open file "//trim(initfile))

       !--------------------------------------------------------------
       ! Read initial positions from init file
       !--------------------------------------------------------------

       if (p_pe==0) WRITE (*,*)
       if (p_pe==0) WRITE (*,*) 'Read initial positions from ', trim(initfile)
       if (p_pe==0) WRITE (*,*)
     
       ! read latitudes
       rcode = nf90_inq_varid (ncid,'LAT',varid)
       IF (rcode /= 0) CALL error_bi('Cannot find variable LAT ',substr)
       rcode = nf90_get_var (ncid,varid,LAT(1:dnparts),start=(/idx(p_pe,1)/), count=(/dnparts/))
       IF (rcode /= 0) CALL error_bi('Cannot read variable LAT ',substr)
       LAT(dnparts+1:dnparts_max) = mdi

       ! read longitude
       rcode = nf90_inq_varid (ncid,'LON',varid)
       IF (rcode /= 0) CALL error_bi('Cannot find variable LON ',substr)
       rcode = nf90_get_var (ncid,varid,LON(1:dnparts),start=(/idx(p_pe,1)/), count=(/dnparts/))
       IF (rcode /= 0) CALL error_bi('Cannot read variable LON ',substr)
       LON(dnparts+1:dnparts_max) = mdi 

       ! read presure
       rcode=nf90_inq_varid(ncid,'PRESS_INIT',varid)
       rcode = nf90_get_var(ncid,varid,LEV(1:dnparts),start=(/idx(p_pe,1)/), count=(/dnparts/))
       IF (rcode /= 0) CALL error_bi('Cannot read variable  '//TRIM(init_vertcoorname),substr)
       LEV(dnparts+1:dnparts_max) = mdi 

       if (lstart) then
        ALLOCATE(PRESSI_INIT(NLON,NLEV+1,NGL))   ! Pa
          ! echam pressure field not available here 
          DO JG=1,NGL
           DO JK=1,NLEV+1
             DO JL=1,NLON
                PRESSI_INIT(JL,JK,JG)= VCT(JK) + VCT(NVCLEV+JK)*APZERO
             END DO
           END DO
          END DO

          CALL clams_locate_parcel_cells(1)

       ELSE
          CALL clams_locate_parcel_cells(2)
       ENDIF
       !
       ! Close init file
       rcode=nf90_close(ncid)    
       ENDIF ! .not. lresume
#endif
  END SUBROUTINE clams_init_tracer
! op_sb_20200226-
!--------------------------------------------------------------------
!ju_ec_20180626+
  SUBROUTINE clams_new_tracer

    !
    ! ju_ec 20180626
    !
    ! This subroutine creates the tracer objects for the CLaMS TRACER set for
    ! parcel data. Tracer object creation is not required for the gridded CLaMS
    ! TRACER set, because these tracer objects are copied to that set in
    ! messy_main_tracer_bi.setup_tracer_set_cl, by the copy_tracer_set command.
    ! 
    ! If you want to tag the tracers with more meta-information, this is the
    ! right routine for that. I haven't done that beyond adding in the units of
    ! the tracers.
    !

  USE messy_main_blather_bi,      ONLY: start_message_bi&
                                       ,end_message_bi&
                                       ,info_bi&
                                       ,error_bi
  USE messy_main_tracer_mem_bi,   ONLY: CLTRSTR, CLGPTRSTR
  USE messy_main_tracer_tools_bi, ONLY: tracer_halt
  USE messy_main_tracer,          ONLY: new_tracer
  USE messy_clams_global,         ONLY: rank, nspec, SPECARR&
                                       ,n_cltr,cl_grid_tracers&
                                       ,clams_grid_verbose

  IMPLICIT NONE

  INTEGER :: status
  INTEGER :: i1, i2
  CHARACTER(LEN=*), PARAMETER :: substr = 'clams_new_tracer'

  CALL start_message_bi( modstr, "Start of clams_new_tracer", substr)

  DO i1 = 1, n_cltr

    CLTR_SPECARR_IND(i1) = -1

    ! Information
    CALL info_bi( "CLTR_NAME = "//cl_grid_tracers(i1)%name, substr )

    ! Find specarr index within SPECARR array
    if (TRIM(cl_grid_tracers(i1)%name).NE.('count')) then

        DO i2 = 1, nspec
            IF (TRIM(SPECARR(i2)%name)==TRIM(cl_grid_tracers(i1)%name))&
                CLTR_SPECARR_IND(i1) = i2
        ENDDO

    else

        CLTR_SPECARR_IND(i1) = 0

    endif

    
    CALL info_bi( "CLTR_SPECARR_IND = "//char(CLTR_SPECARR_IND(i1)), substr)
    if ((rank .EQ. 0) .AND. (CLTR_SPECARR_IND(i1) .EQ. -1)) &
        CALL error_bi("Tracer could not be found in specarr.",substr)

    ! Make the tracer object
    CALL info_bi(&
         "Creating tracer objects for"//cl_grid_tracers(i1)%name&
        ,substr&
        )

    !
    ! ju_ec 20190206
    ! Optional arguments for the new_tracer call are:
    ! unit      (string)
    ! subname   (string)
    ! longname  (string)
    ! medium    (integer)
    ! quantity  (integer)
    ! type      (integer)
    ! ...as well as a few others. See the TRACER manual for details. That can
    ! probably be found by checking the original TRACER paper, which is Joeckel
    ! et al 2008 from ACP, Technical Note: Coupling of chemical processes with
    ! the Modular Earth Submodel System (MESSy) submodel TRACER.
    !

    CALL new_tracer( status, CLTRSTR, cl_grid_tracers(i1)%name, modstr&
        ,idx=CLTR_IDX(i1)&      ! Sets the index of the tracer in the tracer set
        ,unit="mol/mol"&        ! All clams tracer data is mol/mol
        )
    CALL tracer_halt( substr, status )

    !
    ! ju_ec 20190206
    ! One could also use the subroutine set_tracer to set information about the
    ! tracers within the tracer set. This information can be arbitrary, and I
    ! would guess that the information would then be added to the channel output
    ! for the object. One could also access that information via a call to
    ! get_tracer, which could clean up the communication of CLaMS MESSy-TRACER
    ! information. For example, one could place the "CLaMS domain" level
    ! boundaries of each tracer in a string via set_tracer, and then access that
    ! information later in the run, when required.
    !
    
    !
    ! ju_ec 20190206+
    !
    ! I have not finished this code yet.
    !
    ! Check if the tracer is H2O, and if it has a fill channel, and if that
    ! channel has units of kg/kg, in which case the fill channel data must be
    ! multiplied by a molar mass factor before filling the gridded clams data
    ! or the parcel clams data.
    !
    !
    !if (        (TRIM(cl_grid_tracers(i1)%name).EQ.('H2O'))&
    !      .AND. (TRIM(cl_grid_tracers(i1)%fill_ch).NE.('NONE'))) then
    !    
    !    !get channel
    !    !check units
    !    !assign fill factor
    !else
    !    !assign fill factor 1
    !endif
    !CALL set_tracer(status, CLTRSTR, idx, S_fill_factor, fill_factor)
    !CALL tracer_halt( substr, status )
    !CALL set_tracer(status, CLTRSTR, idx, S_fill_channel, fill_factor)
    !CALL tracer_halt( substr, status )
    !CALL set_tracer(status, CLTRSTR, idx, S_fill_object, fill_object)
    !CALL tracer_halt( substr, status )
    !
    ! Actually, I need to compute the fill factor on each step. That will be
    ! based on the molar mass of dry air and the amount of water in each grid
    ! point.
    !
    ! ju_ec 20190206-
    !

  END DO
  
  CALL end_message_bi( modstr, "End of clams_new_tracer.", substr)

  END SUBROUTINE clams_new_tracer
!ju_ec_20180626-
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clams_tracer_update

  !
  ! Author: Edward Charlesworth (IEK-7, Juelich)
  ! Written: July 16, 2018
  !
  ! This subroutine updates the tracer data for the clams parcel
  ! tracer set and calculates the data for the clams gridded
  ! representation.
  !
  ! 
  !

  USE messy_main_tracer_mem_bi, ONLY: xt_c, xt_clgp, CLTRSTR,&
                                      CLGPTRSTR 
  USE messy_clams_global,       ONLY: rank, SPECARR, &
                                      E5_LAT, E5_LON, E5_LEVEL, &
                                      dnparts_max, n_cltr,&
                                      clams_grid_verbose,&
                                      cl_grid_tracers

  USE messy_main_blather_bi,       ONLY: start_message_bi,end_message_bi,info_bi
  USE messy_main_channel,          ONLY: get_channel_object
  USE messy_main_channel_error_bi, ONLY: channel_halt
  USE messy_main_transform_bi,     ONLY: trp_gpdc_gpgl, M_SUM
#ifndef NOMPI
  USE mo_mpi, ONLY: p_all_comm, p_error, MPI_INTEGER, MPI_SUM
#endif

  IMPLICIT NONE

  INTEGER :: status
  CHARACTER(LEN=*), PARAMETER :: substr = 'clams_tracer_update'

  INTEGER :: i1, jt

  INTEGER :: j1,j2,j3
  REAL(DP):: v
  CHARACTER(LEN=20) :: ch, ob
  REAL(DP), DIMENSION(:,:,:), POINTER :: fill => NULL()

  CALL start_message_bi( modstr, "Beginning clams_tracer_update", substr)
  ! Initialize Variables
  pc_d(:,:,:) = 0               ! We must reset these on each call
  pc_g(:,:,:) = 0               !
  clgp_data(:,:,:,:) = 0.0_dp   !
  
! op_sb_20191014+
  ! Update MESSy-TRACER Lagrangian representation
!!$  DO i1 = 1,n_cltr
!!$    jt = CLTR_IDX(i1)
!!$    xt_c(:,1,jt,1) = SPECARR(CLTR_SPECARR_IND(i1))%values
!!$  ENDDO
! op_sb_20191014-

  IF (clams_grid_verbose) CALL info_bi(&
      "Running clams_locate_parcel_cells."&
      ,substr&
      )

  CALL clams_locate_parcel_cells(2) ! Sets the grid-box position for all parcels
                                  ! and counts the number of parcels in each box
                                  ! onto pc_d.

  IF (clams_grid_verbose) CALL info_bi(&
      "Finished clams_locate_parcel_cells."&
      ,substr&
      )

!!$  ! Set up the global grids
!!$  CALL MPI_BARRIER( p_all_comm, p_error ) ! A barrier to ensure that the 
!!$                                          ! allreduce call passes correctly.
!!$  CALL MPI_ALLREDUCE( pc_d, pc_g,& ! Get the global parcels-in-grid-box count
!!$      nlat*nlev*nlon, MPI_INTEGER, MPI_SUM, p_all_comm, p_error )
!!$  CALL MPI_BARRIER( p_all_comm, p_error ) ! A barrier to ensure that the 
!!$                                          ! allreduce call passes correctly.
!!$
  IF (clams_grid_verbose) CALL info_bi(&
      "Finished clams_locate_parcel_cells."&
      ,substr&
      )

  IF (clams_grid_verbose) write(*,*) "Sum of pc_g: ",sum(pc_g)

  IF (clams_grid_verbose) CALL info_bi(&
      "Assigned global parcel count."&
      ,substr&
      )
  IF (clams_grid_verbose) CALL info_bi(&
      "About to call lg2gp."&
      ,substr&
      )
  
!!$  CALL clams_gp2lg(n_cltr, xt_c, pos, status)
!!$  CALL clams_lg2gp(n_cltr, CLTR_IDX, xt_c, clgp_data, pc_g, pos, status )
  CALL clams_gp2lg(n_cltr, xt_c,status)
  CALL clams_lg2gp(n_cltr, CLTR_IDX, xt_c, clgp_data, status )
  IF (clams_grid_verbose) CALL info_bi(&
      "Finished lg2gp."&
      ,substr&
      )

  IF (clams_grid_verbose) &
      write(*,'(A17,E12.5E2,A5,E12.5E2)') &
      "clgp_data max(1) ",MAXVAL(clgp_data(:,:,1,:))&
      ," (2) ",MAXVAL(clgp_data(:,:,2,:))

  gp_d => xt_clgp

  IF (clams_grid_verbose) CALL info_bi(&
      "Transforming from gp_g to gp_d."& 
      ,substr&
      )

  CALL trp_gpdc_gpgl(-1, gp_d, clgp_data, M_SUM )   ! Gather global data

  IF (clams_grid_verbose) &
      write(*,'(A12,E12.5E2,A5,E12.5E2)') &
      "gp_d max(1) ",MAXVAL(gp_d(:,:,1,:))&
      ," (2) ",MAXVAL(gp_d(:,:,2,:))

  do i1 = 1, n_cltr
    ch = cl_grid_tracers(i1)%fill_ch
    ob = cl_grid_tracers(i1)%fill_ob
    if (trim(ch).EQ."NONE") then
      gp_d(:,:,i1,:) = MERGE( gp_d(:,:,i1,:), real(-3.0,dp), gp_d(:,:,i1,:)>0.0 )
    else
      CALL get_channel_object(status, trim(ch), trim(ob), p3=fill)
      CALL channel_halt(substr, status)
      gp_d(:,:,i1,:) = MERGE( gp_d(:,:,i1,:), fill, gp_d(:,:,i1,:)>0.0 )
    end if
  END DO

  IF (clams_grid_verbose) CALL info_bi(&
      "Assigning gp_d to xt_clgp."& 
      ,substr&
      )

  IF (clams_grid_verbose) THEN
    DO i1 = 1, n_cltr

      jt = CLTR_IDX(i1)
      write(*,*) "XT_C      maximum for variable ",TRIM(cl_grid_tracers(i1)%name)&
          ," is ", maxval(xt_c(:,:,jt,:))
      write(*,*) "XT_CLGP   maximum for variable ",TRIM(cl_grid_tracers(i1)%name)&
          ," is ", maxval(xt_clgp(:,:,jt,:))
      write(*,*) "CLGP_DATA maximum for variable ",TRIM(cl_grid_tracers(i1)%name)&
          ," is ", maxval(clgp_data(:,:,jt,:))

    END DO
  ENDIF

  CALL end_message_bi( modstr, "Finished clams_tracer_update", substr)

  END SUBROUTINE clams_tracer_update
!--------------------------------------------------------------------

  SUBROUTINE clams_update_celldist

  !
  ! Author: Edward Charlesworth (IEK-7, Juelich)
  ! Written: July 16, 2018
  !
  ! This subroutine updates the tracer data for the clams parcel
  ! tracer set and calculates the data for the clams gridded
  ! representation.
  !
  ! 
  !

  USE messy_main_tracer_mem_bi, ONLY: xt_c, xt_clgp, CLTRSTR,&
                                      CLGPTRSTR 
  USE messy_clams_global,       ONLY: rank, SPECARR, &
                                      E5_LAT, E5_LON, E5_LEVEL, &
                                      dnparts_max, n_cltr,&
                                      clams_grid_verbose,&
                                      cl_grid_tracers

  USE messy_main_blather_bi,       ONLY: start_message_bi,end_message_bi,info_bi
  USE messy_main_channel,          ONLY: get_channel_object
  USE messy_main_channel_error_bi, ONLY: channel_halt
  USE messy_main_transform_bi,     ONLY: trp_gpdc_gpgl, M_SUM

#ifndef NOMPI
  USE mo_mpi, ONLY: p_all_comm, p_error, MPI_INTEGER, MPI_SUM
#endif

  IMPLICIT NONE

  INTEGER :: status
  CHARACTER(LEN=*), PARAMETER :: substr = 'clams_update_celldist'

  INTEGER :: i1, jt

  INTEGER :: j1,j2,j3
  REAL(DP):: v
  CHARACTER(LEN=20) :: ch, ob
  REAL(DP), DIMENSION(:,:,:), POINTER :: fill => NULL()

  CALL start_message_bi( modstr, "Beginning clams_update_celldist", substr)
  ! Initialize Variables
  pc_d(:,:,:) = 0._dp               ! We must reset these on each call
  pc_g(:,:,:) = 0._dp               
  clgp_data(:,:,:,:) = 0.0_dp   !
  
! op_sb_20191014+
  ! Update MESSy-TRACER Lagrangian representation
!!$  DO i1 = 1,n_cltr
!!$    jt = CLTR_IDX(i1)
!!$    xt_c(:,1,jt,1) = SPECARR(CLTR_SPECARR_IND(i1))%values
!!$  ENDDO
! op_sb_20191014-

  IF (clams_grid_verbose) CALL info_bi(&
      "Running clams_locate_parcel_cells."&
      ,substr&
      )
 
  CALL clams_locate_parcel_cells(2) ! Sets the grid-box position for all parcels
                                  ! and counts the number of parcels in each box
                                  ! onto pc_d.

  IF (clams_grid_verbose) CALL info_bi(&
      "Finished clams_locate_parcel_cells."&
      ,substr&
      )

#ifndef NOMPI
  ! Set up the global grids
  CALL MPI_BARRIER( p_all_comm, p_error ) ! A barrier to ensure that the 
                                          ! allreduce call passes correctly.
  CALL MPI_ALLREDUCE( pc_d, pc_g,& ! Get the global parcels-in-grid-box count
      nlat*nlev*nlon, MPI_INTEGER, MPI_SUM, p_all_comm, p_error )
  CALL MPI_BARRIER( p_all_comm, p_error ) ! A barrier to ensure that the 
                                          ! allreduce call passes correctly.
#endif

  IF (clams_grid_verbose) CALL info_bi(&
      "Finished clams_locate_parcel_cells."&
      ,substr&
      )

  IF (clams_grid_verbose) write(*,*) "Sum of pc_g: ",sum(pc_g)

  IF (clams_grid_verbose) CALL info_bi(&
      "Assigned global parcel count."&
      ,substr&
      )
  IF (clams_grid_verbose) CALL info_bi(&
      "About to call lg2gp."&
      ,substr&
      )
  
!!$  CALL clams_gp2lg(n_cltr, xt_c, pos, status)
!!$  CALL clams_lg2gp(n_cltr, CLTR_IDX, xt_c, clgp_data, pc_g, pos, status )
  CALL clams_gp2lg(n_cltr, xt_c,status)
  CALL clams_lg2gp(n_cltr, CLTR_IDX, xt_c, clgp_data,status )
  IF (clams_grid_verbose) CALL info_bi(&
      "Finished lg2gp."&
      ,substr&
      )

  IF (clams_grid_verbose) &
      write(*,'(A17,E12.5E2,A5,E12.5E2)') &
      "clgp_data max(1) ",MAXVAL(clgp_data(:,:,1,:))&
      ," (2) ",MAXVAL(clgp_data(:,:,2,:))

  gp_d => xt_clgp

  IF (clams_grid_verbose) CALL info_bi(&
      "Transforming from gp_g to gp_d."& 
      ,substr&
      )

  CALL trp_gpdc_gpgl(-1, gp_d, clgp_data, M_SUM )   ! Gather global data

  IF (clams_grid_verbose) &
      write(*,'(A12,E12.5E2,A5,E12.5E2)') &
      "gp_d max(1) ",MAXVAL(gp_d(:,:,1,:))&
      ," (2) ",MAXVAL(gp_d(:,:,2,:))

  do i1 = 1, n_cltr
    ch = cl_grid_tracers(i1)%fill_ch
    ob = cl_grid_tracers(i1)%fill_ob
    if (trim(ch).EQ."NONE") then
      gp_d(:,:,i1,:) = MERGE( gp_d(:,:,i1,:), real(-3.0,dp), gp_d(:,:,i1,:)>0.0 )
    else
      CALL get_channel_object(status, trim(ch), trim(ob), p3=fill)
      CALL channel_halt(substr, status)
      gp_d(:,:,i1,:) = MERGE( gp_d(:,:,i1,:), fill, gp_d(:,:,i1,:)>0.0 )
    end if
  END DO

  IF (clams_grid_verbose) CALL info_bi(&
      "Assinging gp_d to xt_clgp."& 
      ,substr&
      )

  IF (clams_grid_verbose) THEN
    DO i1 = 1, n_cltr

      jt = CLTR_IDX(i1)
      write(*,*) "XT_C      maximum for variable ",TRIM(cl_grid_tracers(i1)%name)&
          ," is ", maxval(xt_c(:,:,jt,:))
      write(*,*) "XT_CLGP   maximum for variable ",TRIM(cl_grid_tracers(i1)%name)&
          ," is ", maxval(xt_clgp(:,:,jt,:))
      write(*,*) "CLGP_DATA maximum for variable ",TRIM(cl_grid_tracers(i1)%name)&
          ," is ", maxval(clgp_data(:,:,jt,:))

    END DO
  ENDIF

  CALL end_message_bi( modstr, "Finished clams_update_celldist", substr)

  END SUBROUTINE clams_update_celldist

!--------------------------------------------------------------------
  SUBROUTINE clams_locate_parcel_cells(flag)
      
      !
      ! Author: Edward Charlesworth (IEK-7, Juelich)
      ! Written: July 17, 2018
      !
      ! This subroutine finds the EMAC grid box that each CLaMS
      ! parcel is located in.
      !
    
 USE messy_main_data_bi,          ONLY: pressi_3d
 USE messy_main_channel,          ONLY: get_channel_object
 USE messy_main_channel_error_bi, ONLY: channel_halt
 USE messy_clams_global,          ONLY: dnparts_max, dnparts 
 USE messy_main_transform_bi,     ONLY: trp_gpdc_gpgl, &
                                        M_SUM   ! op_sb_20190801
 USE messy_main_mpi_bi,           ONLY: p_pe! op_sb_20200304
! op_sb_20200303+
 USE messy_main_grid_def_mem_bi,  ONLY: vct, apzero, nvclev,        &
                                        ngl, nlon, nlev
! op_sb_20200303-

 IMPLICIT NONE

 INTEGER :: status
 CHARACTER(LEN=*), PARAMETER :: substr = 'clams_locate_parcel_cells'
 REAL(dp):: GDLON, GDHLON        ! op_sb_20200407

 INTEGER :: iPar, iLat, iVert, iLon,jl,jk,jg
 ! op_sb_20200303+
 INTEGER, INTENT(IN) :: flag
 ! op_sb_20200303-

 ! Local
 REAL(DP), DIMENSION(:)    , POINTER :: PRESS => NULL()
 INTEGER, PARAMETER :: dim_lon = 1
 INTEGER, PARAMETER :: dim_vrt = 2
 INTEGER, PARAMETER :: dim_lat = 3
 
 ! Longitude grid
 GDLON= 360._dp/REAL(NLON)
 GDHLON= GDLON/2.0_dp

 ! Initialization
 !!$pos(:,:) = -1
 pos(:,:) = -1.      ! op_sb_20190712
 ! Initialize Variables
 pc_d(:,:,:) = 0               ! We must reset these on each call
 pc_g(:,:,:) = 0 

 CALL get_channel_object(status, modstr, 'PRESS', p1=PRESS)
 CALL channel_halt(substr, status)

! op_sb_20200303+
 SELECT CASE (flag)
 CASE(1)   
   ! parcel pressure from init file
   ! GP pressure calculated in init_tracer 
   PRESS(1:dnparts_max) = LEV(1:dnparts_max)   ! in hPa
   CALL trp_gpdc_gpgl(1, pressi_init, pressi_3d_g)
   pressi_3d_g = pressi_3d_g/100.0       
 CASE(2)
   ! Turn pressi_3d into a global field
   ! shape is lon, vert, lat
   CALL trp_gpdc_gpgl(1, pressi_3d, pressi_3d_g) ! Sets pressi_3d_g's shape
   pressi_3d_g = pressi_3d_g/100.0               !sets to hPa, like in CLaMS
 END SELECT
! op_sb_20200303-

 ! Determine bins
 DO iPar = 1, dnparts_max

   if (LAT(iPar)<-90.) CYCLE ! skip undefined parcels - I assume if latitude is 
                             ! undefined then longitude is also undefined, so I
                             ! only check latitude. Same for pressure. Undefined
                             ! here means a very large negative number, but any-
                             ! -thing less than -90 will do here.
   ! ignore parcels above 60 hPa or 300 hPa
   ! if ((PRESS(iPar)< 60.).OR.(PRESS(iPar)>300.)) CYCLE 

   ! Loop over latitudes
   DO iLat = 2,nlat+1 ! We want to loop over the latitude boundaries and we want
                      ! to skip the first one because nothing can have a higher
                      ! latitude than the first (+90)
     if (LAT(iPar)>=E5_LAT_BOUNDS(iLat)) then ! the parcel must be in prev. bin
!!$       pos( dim_lat, iPar ) = iLat-1
       pos(iPar,dim_lat ) = REAL(iLat-1)        ! op_sb_20190712
       EXIT
     endif
   END DO
   
   ! Loop over longitudes
!!$   if (LON(iPar)>E5_LON_BOUNDS(nlon)) then
!!$       pos( dim_lon, iPar ) = 1 ! Checks for parcels east of 0 longitude and in
!!$   else
     DO iLon = 1,nlon
     if (LON(iPar)<=(E5_LON_BOUNDS(iLon)+GDHLON)) then ! parcel must be in prev. bin
!!$       if (LON(iPar)<=E5_LON_BOUNDS(iLon)) then ! parcel must be in prev. bin
!!$         pos( dim_lon, iPar ) = iLon ! prev bin here means the one to the east
         pos(iPar,dim_lon ) = REAL(iLon)  ! op_sb_20190712
         ! do not count parcel position twice (for array spc_g) in case of LON=0.
         if (LON(iPar) == 0._dp) pos(iPar,dim_lon) = nlon ! op_sb_20200407
         EXIT
       endif
     END DO
!!$   endif

  ! Loop over layers
  DO iVert = 1,nlev
    ! some parcel-pressures from initial file are larger then the
    ! GP-surface pressure for parcels in (1:dnparts). This leads to pos=-1
    ! and thus to an error for pc_d ...
    if (iVert == nlev .and. iPar .le. dnparts .and. press(iPar) .gt. &
          pressi_3d_g(INT(pos(iPar,dim_lon)),iVert+1,INT(pos(iPar,dim_lat)))) &
     then
       pos(iPar,dim_vrt) = REAL(iVert)
    endif
    if ( PRESS(iPar)<=&
!!$      pressi_3d_g(pos(dim_lon,iPar),iVert+1,pos(dim_lat,iPar)) ) then
!!$      pos( dim_vrt, iPar) = iVert
! op_sb_20190712
      pressi_3d_g(INT(pos(iPar,dim_lon)),iVert+1,INT(pos(iPar,dim_lat))) ) then
      pos(iPar,dim_vrt) = REAL(iVert)
! op_sb_20190712
      EXIT
    endif
  END DO

  ! ju_ec 20181016
  ! I just learned that some parcels have pressures which are higher than
  ! the maximum pressure in their EMAC column. Therefore these are invalid
  ! and should have their other position indexes set to -1, since, at 
  ! present, these are used to determine if the parcel position was found
  ! and probably all the position indexes should be invalid if one of them
  ! is. Just in the case that somebody tries to check that based on another
  ! position index. So here I check if any of the position indexes are
  ! invalid and then set them all invalid, if that's the case, and then
  ! cycle to the next parcel.
  !
! op_sb_20190712
!!$  if (minval(pos(:,iPar)).LT.1) then
!!$      pos(:,iPar)=-1
  if (minval(pos(iPar,:)) .LT. 1.) then
      pos(iPar,:)=-1.
      CYCLE
  endif

  ! Assign to parcel count
! op_sb_20190712
!!$  pc_d(pos(dim_lon,iPar),pos(dim_vrt,iPar),pos(dim_lat,iPar))&
!!$  = pc_d(pos(dim_lon,iPar),pos(dim_vrt,iPar),pos(dim_lat,iPar)) + 1
  if (pos(iPar,dim_lon) .ge. 1. .and. pos(iPar,dim_vrt) .ge. 1.      &
                                 .and. pos(iPar,dim_lat) .ge. 1.) then ! op_sb_20190802
  pc_d(INT(pos(iPar,dim_lon)),INT(pos(iPar,dim_vrt)),INT(pos(iPar,dim_lat)))&
    = pc_d(INT(pos(iPar,dim_lon)),INT(pos(iPar,dim_vrt)),            &
                                  INT(pos(iPar,dim_lat))) + 1.
  endif ! op_sb_20190802

  ENDDO  !  parcel loop

  CALL trp_gpdc_gpgl(-1, spc_g, pc_d, M_SUM )  ! op_sb_20190802
  CALL trp_gpdc_gpgl(1, spc_g, pc_g )          ! op_sb_20190802

 END SUBROUTINE clams_locate_parcel_cells
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clams_gp2lg(n_cltr,xt_c,status) ! op_sb_20190712
!!$  SUBROUTINE clams_gp2lg(n_cltr,xt_c,pos,status)

    !-----------------------------------------------------------
    !ju_ec_20190213
    !
    ! This subroutine replaces clams parcel data (in SPECARR)
    ! with data from the defined tracer fill channel object, if
    ! one is defined. The replacement only happens to parcels
    ! outside of the clams radiative-feedback coupling region.
    !
    !-----------------------------------------------------------

  USE messy_clams_global,          ONLY: cl_grid_tracers, SPECARR
  USE messy_main_channel,          ONLY: get_channel_object
  USE messy_main_channel_error_bi, ONLY: channel_halt
  USE messy_main_transform_bi,     ONLY: trp_gpdc_gpgl

  IMPLICIT NONE

  !INTRINSIC :: ALLOCATED

  ! I/O
  ! lagrangian field
  REAL(dp), DIMENSION(:,:,:,:), INTENT(INOUT)   :: xt_c
  ! position of parcels
!!$  INTEGER, DIMENSION(:,:), INTENT(IN)           :: pos 
  ! numbers of tracers
  INTEGER, INTENT(IN)                           :: n_cltr
  ! error status
  INTEGER, INTENT(OUT),  OPTIONAL               :: status

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'clams_gp2lg'
  CHARACTER(LEN=20) :: ch, ob
  INTEGER  :: lb, ub
  INTEGER  :: jn, it
  INTEGER  :: i1, i2, i3   ! INDEX IN EXTERNAL GP-FIELDS
  REAL(DP), DIMENSION(:,:,:) , POINTER :: FILL_D => NULL()
  REAL(DP), DIMENSION(:,:,:) , POINTER :: FILL_G => NULL()
  ! LOCAL -- only for conversion of ECHAM5/qm1 channel object to mol/mol
  REAL(DP), PARAMETER :: mmDryAir = 28.9647 ! g/mol
  REAL(DP), PARAMETER :: mmWater  = 18.0153 ! g/mol



  DO it = 1, n_cltr
  IF (cl_grid_tracers(it)%couple) THEN

    lb = cl_grid_tracers(it)%lower  ! lower boundary index
    ub = cl_grid_tracers(it)%upper  ! upper boundary index

    ! Check the fill channel and fill object
    ch = cl_grid_tracers(it)%fill_ch    ! fill channel
    ob = cl_grid_tracers(it)%fill_ob    ! fill channel object

    if (trim(ch).NE."NONE") then

      ! Find and acquire channel object
      CALL get_channel_object(status, trim(ch), trim(ob), p3=fill_d)
      CALL channel_halt(substr, status)

      ! Check if the tracer is water and fill object is qm1 from ECHAM5 and
      ! convert fill_d to mol/mol if that's the case.
      if (     (TRIM(cl_grid_tracers(it)%name).EQ.'H2O')                    &
                .AND.(TRIM(ch).EQ.'ECHAM5')                                 &
                .AND.(TRIM(ob).EQ.'qm1')  )                                 &
         fill_d = fill_d * (                                                &
             mmDryAir / ( 1 + fill_d * ( (mmDryAir - mmWater) / mmWater ) ) &
             ) / mmWater    ! Now fill_d is converted to mol/mol. It may be
      ! prudent to apply this calculation not just to the case of qm1 from
      ! ECHAM5 but also to any channel with kg/kg units (which can be checked
      ! with messy_main_channel/get_attribute), with the final mmWater replaced
      ! by a molar mass for the tracer (that can be stored as an attribute of
      ! the tracer) and the rest of the fill_d on the second line replaced with
      ! the qm1 tracer from ECHAM5 (so that the molar mass of air can be
      ! calculated).

      ! transform the decomposed field (fill_d) to global (fill_g)
      CALL trp_gpdc_gpgl(1, fill_d, fill_g)

!!$      DO jn = 1, size(pos,2)
!!$
!!$        if (pos(1,jn)<0) cycle  ! If parcel is not located, skip it

!!$        i2 = pos(2,jn) ! vertical index of parcel cell

!!$        if ((i2.GE.lb).AND.(i2.LE.ub)) cycle ! if parcel in CLaMS region, skip

!!$        i1 = pos(1,jn) ! longitudinal index of parcel cell
!!$        i3 = pos(3,jn) ! latitudinal index of parcel cell
!  op_sb_20190719+
      DO jn = 1, size(pos,2)

        if (INT(pos(jn,1))<0) cycle  ! If parcel is not located, skip it

        i2 = INT(pos(jn,2)) ! vertical index of parcel cell

        if ((i2.GE.lb).AND.(i2.LE.ub)) cycle ! if parcel in CLaMS region, skip

        i1 = INT(pos(jn,1)) ! longitudinal index of parcel cell
        i3 = INT(pos(jn,3)) ! latitudinal index of parcel cell
! op_sb_20190719-

        ! replace the parcel data with the fill channel object's data
        SPECARR(CLTR_SPECARR_IND(it))%values(jn) = fill_g(i1,i2,i3)
      END DO

    endif

  END IF
  END DO
  
  END SUBROUTINE clams_gp2lg
!--------------------------------------------------------------------


!--------------------------------------------------------------------
  SUBROUTINE clams_lg2gp(n_cltr,memiarr,xt_c,gp,status) ! op_sb_20190712
!!$  SUBROUTINE clams_lg2gp(n_cltr,memiarr,xt_c,gp,pc,status)
  
  IMPLICIT NONE

  !INTRINSIC :: ALLOCATED

  ! I/O
  ! index of tracer in memory
  INTEGER, DIMENSION(:), INTENT(IN)                  :: memiarr 
  ! global grid point field ([mol/mol] or [kg/kg])
  REAL(dp), DIMENSION(:,:,:,:), INTENT(INOUT)        :: gp
  ! lagrangian field
  REAL(dp), DIMENSION(:,:,:,:), INTENT(IN)           :: xt_c
  ! number of tracers
  INTEGER, INTENT(IN)                                :: n_cltr
  ! position of parcels
!!$  INTEGER, DIMENSION(:,:), INTENT(IN)                :: pos
  ! parcel count in each grid box
!!$  INTEGER,  DIMENSION(:,:,:), INTENT(IN)             :: pc
  ! error status
  INTEGER,                    INTENT(OUT),  OPTIONAL :: status

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'clams_lg2gp'
  INTEGER  :: jn, it, lt
  INTEGER  :: i1, i2, i3   ! INDEX IN EXTERNAL GP-FIELDS

  ! INIT
  IF (PRESENT(status)) status = 1   ! ERROR

  lt = n_cltr
  
!!$  DO jn = 1, size(pos,2)
!!$    if (pos(1,jn)<0) cycle  ! If parcel is not located, skip it
!!$    i1 = pos(1,jn) ! lon
!!$    i2 = pos(2,jn) ! vert
!!$    i3 = pos(3,jn) ! lat
! op_sb_20190719+
  DO jn = 1, size(INT(pos),1)
    if (INT(pos(jn,1))<0) cycle  ! If parcel is not located, skip it
    i1 = INT(pos(jn,1)) ! lon
    i2 = INT(pos(jn,2)) ! vert
    i3 = INT(pos(jn,3)) ! lat
! op_sb_20190719-
    DO it = 1, lt
      gp(i1,i2,it,i3)=gp(i1,i2,it,i3)&
          +xt_c(jn,1,memiarr(it),1)/pc_g(i1,i2,i3)   ! op_sb_20190730
!          +xt_c(jn,1,memiarr(it),1)/REAL(pc(i1,i2,i3),dp)   ! op_sb_20190730
    END DO

  END DO
  
  IF (PRESENT(status)) status = 0   ! NO ERROR, if everything worked
  
  END SUBROUTINE clams_lg2gp
!--------------------------------------------------------------------

#endif


!-----------------------------------------------------------------------
#ifdef ECHAM5
!#include "messy_clams_si.inc" ! op_pj_20160713
#include "messy_main_channel_clams.inc"
#endif
!-----------------------------------------------------------------------

#endif
!**********************************************************************
END MODULE messy_clams_si
!**********************************************************************
