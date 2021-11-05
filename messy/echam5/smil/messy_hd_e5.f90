! **********************************************************************
!
!  Hydrological Discharge model
!
!  MESSy- submodel interface for HD
!
!  AUTHOR:  Pozzer Andrea, MPICH, October 2007
!           The original code is based on the
!           HD model ***** Version 1.0 - Oktober 1999
!           Programmed and developed by Stefan Hagemann, MPIM
!           For specific routine contribution check the code!
! mz_bk_20120730: Update: "Speed-up" (ECHAM6/MPIESM-1.0.00)
! mz_bk_20120923: Disabled "new" glacier calving, using "old" glacier_to_ocean
!                 WARNING: Several changes have to be implemented to use
!                          "new" version!!
! mz_bk_20120928: Switch "lnew_glac" to switch to "new" ECHAM6 glacier calving
!                 .true. = ECHAM6 version ("direct" treatment)
!                 .false.= ECHAM5 version (glacier_to_ocean) (default)
!
! **********************************************************************
MODULE messy_hd_e5

  ! MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_pe, p_bcast,   & 
                                      finish, message
  USE messy_main_channel,       ONLY: REPR_UNDEF, DIMID_UNDEF
  USE messy_main_timer_bi,      ONLY: p_bcast_event,timer_event_init
  USE messy_main_timer_event,   ONLY: io_time_event, TRIG_FIRST      &
                                    , TIME_INC_DAYS, time_event
  USE messy_main_data_bi,       ONLY: lcouple
  USE mo_constants,             ONLY: rhoh2o
  USE messy_hd

  IMPLICIT NONE 
  CHARACTER(512) :: message_text = ''

  ! NEW CHANNEL REPRESENTATION
  INTEGER :: GP_2D_HD     = REPR_UNDEF
  INTEGER :: GP_3D_HD     = REPR_UNDEF
  INTEGER :: DIMID_LON_HD   = DIMID_UNDEF
  INTEGER :: DIMID_LAT_HD   = DIMID_UNDEF
  INTEGER :: DIMID_LEV_HD   = DIMID_UNDEF


  PRIVATE

  ! TIME MANAGER
  LOGICAL, PUBLIC, SAVE ::  l_trig_hd = .TRUE.  
  LOGICAL, PUBLIC, SAVE ::  l_trig_hd_glacier = .TRUE.  
  TYPE(io_time_event), PUBLIC, SAVE :: trig_hd = &
       io_time_event (1,'days',TRIG_FIRST,0)
  TYPE(io_time_event), PUBLIC, SAVE :: trig_hd_glacier = & 
       io_time_event (1,'days',TRIG_FIRST,0)
  TYPE(time_event),    PUBLIC, SAVE :: ev_trig_hd
  TYPE(time_event),    PUBLIC, SAVE :: ev_trig_hd_glacier
  ! op_bk_20161020+
!!$ INTEGER, SAVE :: hd_interval 
  REAL(DP), SAVE :: hd_interval 
  ! op_bk_20161020-
  INTEGER, SAVE :: hd_glacier_interval 

  NAMELIST /CPL/ trig_hd, trig_hd_glacier

  ! mz_bk_20120730+
  TYPE cart_idx_2d
     INTEGER :: ilon, ilat
  END TYPE cart_idx_2d

  TYPE cart_coord_2d
     REAL(dp) :: lon, lat
  END TYPE cart_coord_2d

  TYPE cart_xidx_2d
     INTEGER :: ilon, ilat, extlen
     REAL(dp), ALLOCATABLE :: amod(:), akdiv(:)
  END TYPE cart_xidx_2d
  ! mz_bk_20120730-

  ! HD model grid dimensions
  INTEGER, PARAMETER :: nl = 720    ! number of longitudes
  INTEGER, PARAMETER :: nb = 360    ! number of latitudes
  ! mz_bk_20120730+
  REAL(dp), PARAMETER :: fullcirc = 360.0_dp
  REAL(dp), PARAMETER :: halfcirc = fullcirc * 0.5_dp
  REAL(dp), PARAMETER :: hd_scal_lon = fullcirc/nl
  REAL(dp), PARAMETER :: hd_scal_lat = 0.5_dp*fullcirc/nb
  !
  ! mm: Sub (computational) time steps per day for Riverflow = 4
  INTEGER, PARAMETER :: mm = 4
  !
  ! ndd = If no land point is found as direct neighbour,
  !       it is searched in NWSE direction until the maximum distance of
  !       NDD Boxes is reached.
  INTEGER, PARAMETER :: ndd = 3
  !
  ! flow types, used in kasglob
  INTEGER, PARAMETER :: overlandflow = 1, riverflow = 2
  ! mz_bk_20120730-

  ! North-west corner of gridbox(1,1) and resolution

  REAL(dp), PARAMETER :: florg = -180.0_dp
  REAL(dp), PARAMETER :: fborg =   90.0_dp
  REAL(dp), PARAMETER :: fscal =    0.5_dp

  ! mz_bk_20120730+
  ! corresponding coordinates on echam ocean grid
  REAL(dp), PARAMETER :: oclorg = 0.0_dp, ocborg = 90.0_dp
  ! mz_bk_20120730-

  INTEGER, PARAMETER :: nmemrf = 5

  REAL(dp), POINTER :: alf_k(:,:)    => NULL()! retention constant k, overflow
  REAL(dp), POINTER :: alf_n(:,:)    => NULL()! number of reservoirs n , overflow
  REAL(dp), POINTER :: arf_k(:,:)    => NULL()! retention constant k, riverflow
  REAL(dp), POINTER :: arf_n(:,:)    => NULL()! number of reservoirs  n , riverflow
  REAL(dp), POINTER :: agf_k(:,:)    => NULL()! retention constant k, baseflow
  ! mz_bk_20120730+
  ! REAL(dp), POINTER :: fdir(:,:)     => NULL()! river direction
  INTEGER, POINTER :: fdir(:,:)     ! river direction
  ! REAL(dp), POINTER :: flag(:,:)     => NULL()! land mask
  REAL(sp), POINTER :: hd_lsm(:,:)     ! land mask
  ! mz_bk_20120730-
  ! mz_bk_20120730+
  ! REAL(dp), POINTER :: friv(:,:)     => NULL()
  ! REAL(dp), POINTER :: finp(:,:)     => NULL()
  ! REAL(dp), POINTER :: fdata(:,:)    => NULL()
  ! mz_bk_20120730-
  REAL(dp), POINTER :: area(:)       => NULL()
  REAL(dp), POINTER :: finfl(:,:)    => NULL()! Inflow data
  REAL(dp), POINTER :: fgmem(:,:)    => NULL()! intermediate linear baseflow reservoir
  ! mz_bk_20120730+
  TYPE(cart_idx_2d), ALLOCATABLE :: oclook_cache(:,:), intpol_mapping(:,:)
  TYPE(cart_xidx_2d) , ALLOCATABLE :: arf_n_kas(:), alf_n_kas(:)
  ! mz_bk_20120730-
  REAL(dp), POINTER :: frfmem(:,:,:) => NULL()! intermediate reservoirs, inflow cascade
  REAL(dp), POINTER :: flfmem(:,:,:) => NULL()! intermediate reservoir, linear overflow

  REAL(dp), ALLOCATABLE, TARGET :: gl_aros(:,:)
  REAL(dp), ALLOCATABLE, TARGET :: gl_adrain(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_disch(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_zcalv(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_slm(:,:)  
  REAL(dp), ALLOCATABLE, TARGET :: gl_slf(:,:)  
  REAL(dp), ALLOCATABLE, TARGET :: gl_alake(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_awfre(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_aifre(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_awhea(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_aicon(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_apmecal(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_glac(:,:) 
  !mz_ap_20090306+
  REAL(dp), ALLOCATABLE, TARGET :: zpmeb(:,:) 
  !
  ! POINTERS FOR EXTERNAL CHANNEL OBJECTS (replacing USE from messy_main_data_bi)
  REAL(dp), DIMENSION(:,:),   POINTER :: ahfice     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: awhea     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: awfre     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: aifre     => NULL()
  !
  ! see mo_memory_g3b
  !
  !  variables for coupling with HD-model and calving model only
  !
  REAL(dp), POINTER :: acc_aros(:,:)
  REAL(dp), POINTER :: acc_adrain(:,:)
  REAL(dp), POINTER :: aros(:,:)
  REAL(dp), POINTER :: adrain(:,:)
  REAL(dp), POINTER :: apmecal(:,:)
  REAL(dp), POINTER :: disch(:,:)
  REAL(dp), POINTER :: disch_m3s(:,:)
  REAL(dp), POINTER :: zcalv(:,:)   
  REAL(dp), POINTER :: awfre_acc(:,:)
  REAL(dp), POINTER :: aifre_acc(:,:)

  ! water flux correction if coupled model
  REAL(dp), POINTER :: zcorr
  REAL(dp), POINTER :: apmebco(:,:) => NULL()
  REAL(dp), ALLOCATABLE :: qtold(:,:) 
  REAL(dp), ALLOCATABLE :: qtnew(:,:) 

  
  INTRINSIC TRIM

  PUBLIC :: hd_initialize
  PUBLIC :: hd_init_memory
  PUBLIC :: hd_init_coupling
  PUBLIC :: hd_global_end
  PUBLIC :: hd_free_memory


CONTAINS

  ! ---------------------------------------------------------------------------

  SUBROUTINE hd_initialize

    USE messy_main_constants_mem,      ONLY: api => pi, a=>radius_earth
    USE messy_main_grid_def_mem_bi,    ONLY: nlon, ngl 
    USE messy_main_channel_error_bi,   ONLY: channel_halt
    USE messy_main_channel_bi,         ONLY: DC_BC
    USE messy_main_channel_dimensions, ONLY: new_dimension,              &
                                             add_dimension_variable_att, &
                                             add_dimension_variable
    USE messy_main_channel_repr,       ONLY: new_representation
    USE messy_main_tools,              ONLY: find_next_free_unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='hd_initialize'         


    ! HD original
    REAL(dp) :: ra, rb, rh, rd

    INTEGER                 :: status
    INTEGER                 :: iou    ! I/O unit

    ! new rapresentation related
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: array

    INTEGER :: i

    CALL start_message_bi(modstr, 'INITIALISATION', substr)

    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CORE ROUTINE:
       CALL hd_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL finish(substr)
    END IF

    ! BROADCAST RESULTS
    CALL p_bcast(inithd_files_path, p_io)
    CALL p_bcast(lhd_que, p_io)
    CALL p_bcast(lwater_corr, p_io)
    CALL p_bcast(lnew_glac, p_io) ! mz_bk_20120928

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL hd_read_nml_cpl(status, iou)
       IF (status /= 0) CALL finish(substr)
    END IF

    CALL p_bcast_event (trig_hd, p_io)
    CALL p_bcast_event (trig_hd_glacier, p_io)

    ! Initialize memory for the HD Model
    ! corresponds to offline routine 'hdini.f' by S. Hagemann

    IF (p_pe == p_io) THEN

      ALLOCATE (alf_k(nl,nb))          ; alf_k(:,:)    = 0.0_dp
      ALLOCATE (alf_n(nl,nb))          ; alf_n(:,:)    = 0.0_dp
      ALLOCATE (arf_k(nl,nb))          ; arf_k(:,:)    = 0.0_dp
      ALLOCATE (arf_n(nl,nb))          ; arf_n(:,:)    = 0.0_dp
      ALLOCATE (agf_k(nl,nb))          ; agf_k(:,:)    = 0.0_dp
      ! mz_bk_20120730+
      ! ALLOCATE (fdir(nl,nb))           ; fdir(:,:)     = 0.0_dp
      ALLOCATE (fdir(nl,nb))           ; fdir(:,:)     = 0
      ! ALLOCATE (flag(nl,nb))           ; flag(:,:)     = 0.0_dp
      ALLOCATE (hd_lsm(nl,nb))         ; hd_lsm(:,:)   = 0.0_sp
      ! mz_bk_20120730-
      ! mz_bk_20120730+
      ! ALLOCATE (friv(nl,nb))           ; friv(:,:)     = 0.0_dp
      ! ALLOCATE (finp(nl,nb))           ; finp(:,:)     = 0.0_dp
      ! ALLOCATE (fdata(nl,nb))          ; fdata(:,:)    = 0.0_dp
      ! mz_bk_20120730-
      ALLOCATE (area(nb))              ; area(:)       = 0.0_dp
! output (needed for rerun)
!      ALLOCATE (finfl(nl,nb))          ; finfl(:,:)    = 0.0_dp
!      ALLOCATE (fgmem(nl,nb))          ; fgmem(:,:)    = 0.0_dp
!      ALLOCATE (frfmem(nl,nb,nmemrf))  ; frfmem(:,:,:) = 0.0_dp
!      ALLOCATE (flfmem(nl,nb,1))       ; flfmem(:,:,:) = 0.0_dp
      ! mz_bk_20120730+
      ALLOCATE (oclook_cache(nl,nb))   ; oclook_cache  = cart_idx_2d(-1, -1)
      ALLOCATE (intpol_mapping(nl, nb)) ; intpol_mapping = cart_idx_2d(-1, -1)
      ! mz_bk_20120730-
    END IF
    ALLOCATE (gl_aros(nlon,ngl))    ; gl_aros(:,:)    = 0.0_dp
    ALLOCATE (gl_adrain(nlon,ngl))  ; gl_adrain(:,:)  = 0.0_dp  
    ALLOCATE (gl_disch(nlon,ngl))   ; gl_disch(:,:)   = 0.0_dp
    ALLOCATE (gl_zcalv(nlon,ngl))   ; gl_zcalv(:,:)   = 0.0_dp
    ALLOCATE (gl_slm(nlon,ngl))     ; gl_slm(:,:)     = 0.0_dp
    ALLOCATE (gl_slf(nlon,ngl))     ; gl_slf(:,:)     = 0.0_dp
    ALLOCATE (gl_alake(nlon,ngl))   ; gl_alake(:,:)   = 0.0_dp
    IF (lcouple) THEN
      ALLOCATE (gl_awfre(nlon,ngl))   ; gl_awfre(:,:)   = 0.0_dp
      ALLOCATE (gl_aifre(nlon,ngl))   ; gl_aifre(:,:)   = 0.0_dp
      ALLOCATE (gl_awhea(nlon,ngl))   ; gl_awhea(:,:)   = 0.0_dp
    ENDIF
    ALLOCATE (gl_aicon(nlon,ngl))   ; gl_aicon(:,:)   = 0.0_dp
    ALLOCATE (gl_apmecal(nlon,ngl)) ; gl_apmecal(:,:) = 0.0_dp 
    ALLOCATE (gl_glac(nlon,ngl))    ; gl_glac(:,:)    = 0.0_dp
    !mz_ap_20090306+
    ALLOCATE (zpmeb(nlon,ngl));  zpmeb(:,:)   = 0.0_dp

    IF (p_pe == p_io) THEN

      ! Read parameter fields and restart file for the HD Model
      
      !CALL read_hydrology (done in init memory)

      ! setup area in m^2 of the HD model internal grid
    
      ra = 2.0_dp*api*a*a/REAL(nl,dp)
      rb = api/REAL(nb,dp)
      rd = 0.5_dp*api

      DO i = 1, nb
        rh = SIN(-rd+(i-1)*rb)-SIN(-rd+i*rb)
        area(i) = ABS(rh)*ra 
      END DO

    END IF
    ! mz_bk_20120730+
!    CALL hydrology_slm_invariants
    ! mz_bk_20120730-


!-----------------------------------------------------------------------
!                      NEW REPRESENTATION
!-----------------------------------------------------------------------

    ! (1b) NEW DIMENSIONS
    !-----------------  LON -----------------------------------

    ALLOCATE(array(nl))
    DO i=1, nl
       array(i) = -180.0_dp + 180.0_dp/nl + (REAL(i,DP)-1)*360.0_dp/nl
    END DO

    CALL new_dimension(status, DIMID_LON_HD, 'hd_lon', nl)
    CALL channel_halt(substr, status)

    CALL add_dimension_variable(status, 'hd_lon', 'hd_lon',array )
    CALL channel_halt(substr, status)

    CALL add_dimension_variable_att(status, 'hd_lon', 'hd_lon', &
         'long_name', c='hd model longitude')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'hd_lon', 'hd_lon', &
         'units', c='degree east')
    CALL channel_halt(substr, status)

    DEALLOCATE(array)

    !-----------------  LAT -----------------------------------

    ALLOCATE(array(nb))
    DO i=1, nb
       array(i) = 90.0_dp - 90.0_dp/nb - (REAL(i,DP)-1)*180.0_dp/nb
    END DO

    CALL new_dimension(status, DIMID_LAT_HD, 'hd_lat', nb)
    CALL channel_halt(substr, status)

    CALL add_dimension_variable(status, 'hd_lat', 'hd_lat',array )
    CALL channel_halt(substr, status)

    CALL add_dimension_variable_att(status, 'hd_lat', 'hd_lat', &
         'long_name', c='hd model latitude')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'hd_lat', 'hd_lat', &
         'units', c='degree north')
    CALL channel_halt(substr, status)

    DEALLOCATE(array)

    !-----------------  LEV -----------------------------------

    ALLOCATE(array(nmemrf))
    DO i=1, nmemrf
       array(i) = REAL(i,DP)
    END DO

    CALL new_dimension(status, DIMID_LEV_HD, 'hd_lev', nmemrf)
    CALL channel_halt(substr, status)

    CALL add_dimension_variable(status, 'hd_lev', 'hd_lev',array )
    CALL channel_halt(substr, status)

    CALL add_dimension_variable_att(status, 'hd_lev', 'hd_lev', &
         'long_name', c='hd model level')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'hd_lev', 'hd_lev', &
         'units', c='')
    CALL channel_halt(substr, status)

    DEALLOCATE(array)

    !-------------------------------------------------

    ! (1c) NEW REPRESENTATION
    CALL new_representation(status, GP_2D_HD, 'GP_2D_HD' &
         , rank = 2, link = 'xx--', dctype = DC_BC       &
         , dimension_ids = (/ DIMID_LON_HD, DIMID_LAT_HD /)    &
         , ldimlen       = (/ nl , nb /)                       &
         , output_order  = (/ 1,2 /)                           &
         , axis = 'XY--'                                       &
         )
    CALL channel_halt(substr, status)
    CALL new_representation(status, GP_3D_HD, 'GP_3D_HD' &
         , rank = 3, link = 'xxx-', dctype = DC_BC       &
         , dimension_ids = (/ DIMID_LON_HD, DIMID_LAT_HD,  DIMID_LEV_HD /)    &
         , ldimlen       = (/ nl , nb, nmemrf /)                       &
         , output_order  = (/ 1,2,3 /)                           &
         , axis = 'XYZ-'                                       &
         )
    CALL channel_halt(substr, status)


    CALL end_message_bi(modstr, 'INITIALISATION', substr)

!-----------------------------------------------------------------------
!                          EVENT INITIALIZATION
!-----------------------------------------------------------------------

  CALL timer_event_init(ev_trig_hd, trig_hd, &
       'hd computation', 'present')
  CALL timer_event_init(ev_trig_hd_glacier, trig_hd_glacier, &
       'hd computation', 'present')

  IF (p_pe==p_io) WRITE (*,*) 'trig_hd: '         ,trig_hd
  IF (p_pe==p_io) WRITE (*,*) 'trig_hd_glacier: ' ,trig_hd_glacier

  END SUBROUTINE hd_initialize 

  ! ---------------------------------------------------------------------------

  SUBROUTINE hd_init_memory

    ! ***** Version 1.0 - Oktober 1999
    !            Programmed and developed by Stefan Hagemann, MPI
    !
    !            Remark: Input data of Runoff and Drainage should have the
    !                       unit m/s.
    !
    ! ***** Version 1.1 - January 2001
    !       ECHAM5-Version
    !
    ! S.Legutke MPI M&D, Jan 2002, deallocate variables at end of
    !                              rerun cycle
    !
    ! **** MESSy version
    !      
    ! A.Pozzer 2007
    !
    ! Reads parameter fields and restart file for the HD Model
    !
    !       area = array of gridbox areas, Unit = [m^2]
    !    
    !    **** file names
    !    
    !      yodnres = Restart file with reservoir cascade arrays,...
    !      yodnpar = Parameter file with Landmask, RDF, ...

    ! ****** list of variables for rerun
    !
    !   lures = Logical Unit of binar y restart file yodnres
    ! yodnres = Restart file with reservoir cascade arrays,...
    !   nstep = Actual timestep for writing the restart file
    !    ique = Log-Output switch ( 0 = No Log-output to STDOUT)
    !
    !  frfmem(nl, nb, nmemrf) = Intermediate content of reservoir cascade
    !                           for the inflows per Gridbox (=5)
    !
    !  flfmem(nl, nb) = Intermediate content of linear reservoir for
    !                           Overland Flow
    !
    !   fgmem = Array of linear baseflow reservoir (Intermediate content)
    !           At Initialization it has the unit [m^3/s] 
    !               (daily time step inherently implemented)
    !   finfl = Inflow data array for each gridbox for time step nstep

    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: DC_BC, GP_2D_HORIZONTAL, SCALAR


    ! LOCAL
    CHARACTER(len=*), PARAMETER  :: substr = 'hd_init_memory'

    ! CHANNEL MANAGEMENT
    INTEGER                 :: status
    INTEGER                 :: iou    ! I/O unit
    INTRINSIC TRIM


    CALL start_message_bi(modstr, 'MEMORY INITIALIZATION', substr)



! NEW CHANNEL (HD)

    CALL new_channel(status, modstr, reprid=GP_2D_HD, lrestreq=.TRUE.)
    CALL channel_halt(substr,status)

! NEW CHANNEL OBJECT

!--------------------   Linear overlandflow reservoir

    CALL new_channel_object(status, modstr, 'FLFMEM' &
         , p3 = flfmem, reprid = GP_3D_HD)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'FLFMEM' &
                , 'long_name', c=' Linear overlandflow reservoir' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'FLFMEM' &
                , 'units', c='m**3')
    CALL channel_halt(substr, status)

!--------------------    Inflow reservoir cascade

    CALL new_channel_object(status, modstr, 'FRFMEM' &
         , p3 = frfmem, reprid= GP_3D_HD)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'FRFMEM' &
                , 'long_name', c=' Inflow reservoir cascade' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'FRFMEM' &
                , 'units', c='m**3')
    CALL channel_halt(substr, status)

!--------------------    Linear baseflow reservoir

    CALL new_channel_object(status, modstr, 'FGMEM' &
         , p2 = fgmem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'FGMEM' &
                , 'long_name', c=' Linear baseflow reservoir' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'FGMEM' &
                , 'units', c='m**3')
    CALL channel_halt(substr, status)

!--------------------   Inflow for each gridbox 

    CALL new_channel_object(status, modstr, 'FINFL' &
         , p2 = finfl)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'FINFL' &
                , 'long_name', c=' Inflow for each gridbox' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'FINFL' &
                , 'units', c='m**3')
    CALL channel_halt(substr, status)

!--------------------   accumulated atmospheric runoffs 

    CALL new_channel_object(status, modstr, 'acc_aros' &
         , p2 = acc_aros, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'acc_aros' &
                , 'long_name', c='atmospheric runoff accumulated' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'acc_aros' &
                , 'units', c='m')
    CALL channel_halt(substr, status)

!--------------------   accumulated atmospheric drainage

    CALL new_channel_object(status, modstr, 'acc_adrain' &
         , p2 = acc_adrain, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'acc_adrain' &
                , 'long_name', c='atmospheric drainage accumulated' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'acc_adrain' &
                , 'units', c='m')
    CALL channel_halt(substr, status)

!--------------------    atmospheric runoffs 

    CALL new_channel_object(status, modstr, 'aros' &
         , p2 = aros, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'aros' &
                , 'long_name', c='atmospheric runoff' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aros' &
                , 'units', c='m/s')
    CALL channel_halt(substr, status)

!--------------------    atmospheric drainage

    CALL new_channel_object(status, modstr, 'adrain' &
         , p2 = adrain, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'adrain' &
                , 'long_name', c='atmospheric drainage ' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'adrain' &
                , 'units', c='m/s')
    CALL channel_halt(substr, status)

!------------------- atmospheric runoffs 

    CALL new_channel_object(status, modstr, 'apmecal' &
         , p2 = apmecal, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'apmecal' &
                , 'long_name', c='(p - e) at glacier points (accumulated) ' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'apmecal' &
                , 'units', c='m/s')
    CALL channel_halt(substr, status)

    IF (lcouple) THEN
!------------------- accumulated awfre

    CALL new_channel_object(status, modstr, 'awfre_acc' &
         , p2 = awfre_acc, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'awfre_acc' &
                , 'long_name', c='accumulated_awfre' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'awfre_acc' &
                , 'units', c='m')
    CALL channel_halt(substr, status)

!------------------- accumulated awfre

    CALL new_channel_object(status, modstr, 'aifre_acc' &
         , p2 = aifre_acc, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'aifre_acc' &
                , 'long_name', c='accumulated_aifre' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aifre_acc' &
                , 'units', c='m')
    CALL channel_halt(substr, status)

    ENDIF
!------------------- atmospheric total discharge 

    CALL new_channel_object(status, modstr, 'disch' &
         , p2 = disch, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'disch' &
                , 'long_name', c='atmospheric total discharge' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'disch' &
                , 'units', c='m/s')
    CALL channel_halt(substr, status)

!------------------- discharge in m^3s

    CALL new_channel_object(status, modstr, 'disch_m3s' &
         , p2 = disch_m3s, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'disch_m3s' &
                , 'long_name', c='discharge in m3s (without ice melting)' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'disch_m3s' &
                , 'units', c='m3/s')
    CALL channel_halt(substr, status)

!------------------- glacier discharge 

    CALL new_channel_object(status, modstr, 'zcalv' &
         , p2 = zcalv, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'zcalv' &
                , 'long_name', c='glacier discharge' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zcalv' &
                , 'units', c='m3/s')
    CALL channel_halt(substr, status)

!------------------- awfre correction
 
    IF (lwater_corr.and.lcouple) THEN

      CALL new_channel_object(status, modstr, 'zcorr' &
           , p0 = zcorr, reprid=SCALAR)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr , 'zcorr' &
                  , 'long_name', c='water flux into ocean correction' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'zcorr' &
                  , 'units', c='m/s')
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr,  'apmebco', &
           p2=apmebco, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'apmebco', &
           'long_name', c='vert.integr.tendencies of water for correction in hd')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'apmebco', 'units', c='m/s')
      CALL channel_halt(substr, status)

    ENDIF

!------------------------------------------------------------
!------------------------------------------------------------


    IF (p_pe == p_io) THEN

      CALL hd_read_start 

    END IF ! p_pe = p_io

!HB 2012-09-03
    CALL hydrology_slm_invariants
!HB 2012-09-03

    CALL end_message_bi(modstr, 'MEMORY INITIALIZATION', substr)


  END SUBROUTINE hd_init_memory
  ! ---------------------------------------------------------------------

  SUBROUTINE hd_init_coupling

    ! ECHAM5/MESSy
    ! MESSy
    USE messy_main_channel,         ONLY: get_channel_object              & 
                                        , get_channel_object_info         
                                          
    IMPLICIT NONE

    INTRINSIC :: TRIM

    CHARACTER(LEN=*), PARAMETER :: substr = 'hd_init_coupling'
    INTEGER                     :: status

    IF (lcouple) THEN
      CALL get_channel_object(status          &
           , oname='ahfice', cname='ECHAM5'  & 
           ,p2=ahfice)
      IF (status /= 0) THEN 
         WRITE(*,*) "NO channel ahfice from ECHAM5"
         CALL finish(substr)
      ENDIF
      CALL get_channel_object(status          &
           , oname='awhea', cname='a2o'  & 
           ,p2=awhea)
      IF (status /= 0) THEN 
         WRITE(*,*) "NO channel awhea from a2o"
         CALL finish(substr)
      ENDIF
      CALL get_channel_object(status          &
           , oname='awfre', cname='a2o'  & 
           ,p2=awfre)
      IF (status /= 0) THEN 
         WRITE(*,*) "NO channel awfre from a2o"
         CALL finish(substr)
      ENDIF
      CALL get_channel_object(status          &
           , oname='aifre', cname='a2o'  & 
           ,p2=aifre)
      IF (status /= 0) THEN 
         WRITE(*,*) "NO channel aifre from a2o"
         CALL finish(substr)
      ENDIF
    ENDIF

  END SUBROUTINE hd_init_coupling

  ! ---------------------------------------------------------------------------

  SUBROUTINE hd_global_end
  !In ECHAM5 this is called one per day!

    USE messy_main_grid_def_mem_bi, ONLY: nlon, ngl,nn &
                                        , nproma, npromz, kproma, ngpblks  &
                                        , nlev

    USE messy_main_grid_def_bi,     ONLY: gridarea, philat,gboxarea_2d
    USE messy_main_data_bi,         ONLY: slm,slf, alake,glac,          &
                                          aros_2d,adrain_2d,apmecal_2d, &
                                          disch_2d
! for precipitation correction (moved from physc)
    USE messy_main_data_bi,            ONLY: qm1, xlm1, xim1,            &
                                             pressi_3d, rsfl_2d, ssfl_2d,      &
                                             rsfc_2d, ssfc_2d, qflux
     
    USE messy_main_constants_mem,      ONLY: alf, g
    USE messy_main_mpi_bi,             ONLY: scatter_gp, dcg, gather_field
    !for water correction
    USE mo_gaussgrid,                  ONLY: gl_budw
    USE messy_main_timer,              ONLY: delta_time         
    USE messy_main_constants_mem,      ONLY: rho_H2O
    USE messy_main_transform_bi,       ONLY: trp_gpdc_gpgl
    USE messy_main_timer_bi,           ONLY: event_state                   &
                                        , get_time_step => timer_get_time_step
    USE messy_main_timer,              ONLY: current_date, previous_date &
                                           , time_step_len,lstart
    USE messy_main_tools,              ONLY: find_next_free_unit
    

    ! HD Model - Constants and Switches
    ! Also the model grid characteristics of the global ocean grid are added
    !
    ! ibase : Baseflow ON (1) or OUT (0)

    INTEGER, PARAMETER :: ibase = 1

    ! mz_bk_20120730+
    ! iocean (INTEGER) -> locean (LOGICAL)
    ! !
    ! ! iocean: Closure of Water budget for Ocean coupling: 0 = OUT, 1 = ON

    ! INTEGER, PARAMETER :: iocean = 1
    ! locean: Closure of Water budget for Ocean coupling
    LOGICAL, PARAMETER :: locean = .TRUE.
    ! mz_bk_20120730-

    ! mz_bk_20120730+
    ! Moved to definition of module constants
    ! !
    ! ! mm: Sub (computational) time steps per day for Riverflow = 4

    ! INTEGER, PARAMETER :: mm = 4
    ! mz_bk_20120730-

    !
    ! **** Global/Regional Discharge Simulation as Subroutine for ECHAM5
    !
    !
    ! ***** Version 1.0 - November 1999
    !   Programmed and Developed by Stefan Hagemann, MPI
    !   Program code is based on Offline-Version of the HD model
    !   which is also regionally applicable. (regsim.f)
    !
    !   Anmerkung: Input data of Runoff und Drainage should have the unit m/s.
    !
    ! **** Remarks: Changes with regard to offline version
    !   Runoff array is now passed into routine instead of reading it
    !   in kasglob via echread. In echread, now named hdech and called 
    !   before kasglob, only the transformation of the runoff array 
    !   to the resolution of 0.5 degree is done if necessary.
    !   Parameter/Variables luinp, area are deleted from kasglob.
    !
    !   Since the input array to be transformed is only passed to echread
    !   but not read in ECHREAD itself, in ECHREAD only 
    !   the transformation of the input array to 0.5 degree is done.
    !   Therefor, new calling parameter and routine names are set/given.
    !   old: CALL echread(luinp, ihead, tocode, istep, ique)
    !   new: CALL hdech(code_t42, tocode_0.5grad, ique)
    !
    !
    ! ***** River Direction File (RDF) format:
    !
    !                    7  8  9
    !                     \ | /
    !                      \|/
    !                    4--5--6
    !                      /|\
    !                     / | \
    !                    1  2  3
    !
    !       Remark: Direction 5 = Discharge Trap
    !               Direction -1 = Ocean Point
    !
    !
    ! ****** List of variables
    !
    !  ibase = Baseflow ON (1) or OUT (0)
    ! iocean = Closure of Water budget for ocean coupling: 0=Out, 1=On
    !
    ! isolog = Logfile output into Iso file (ASCII file , two columns)
    !      0 = no, 1 = Bothnian Bay/Sea, 2 = Torneaelven, 3 = Global...
    !      4 = St.Lawrence, 5 = Paraguay 6 = Odra
    !lhd_que = Log-Output switch  ( .FALSE. = No Log-Output to STDOUT)
    !
    !  istep = Chosen  time step for Reading of Input
    !     nl = Number of Longitudes
    !     nb = Number of Latitudes
    !
    !     mm = Computation of riverflow in mm internal time steps per day
    !
    ! **** Global Arrays:
    !
    !  finp = local input data array for time step istep
    ! fdata = local output data array for time step istep
    ! finfl = Inflow data array for each gridbox for time step istep
    !  fdir = River direction array
    !hd_lsm = Land mask array
    ! alf_k = Array of retention constants k  - Overland flow [day]
    ! alf_n = Array of number of reservoirs n - Overland flow
    ! arf_k = Array of retention constants k  - Riverflow [day]
    ! arf_n = Array of number of reservoirs n - Riverflow
    ! agf_k = Array of retention constants k  - Baseflow [day]
    !
    !  frfmem(nl, nb, nmemrf) = Intermediate array of reservoir cascade for
    !                           the inflows per Gridbox (new: = nmemrf = 5)
    !
    !  flfmem(nl, nb, nmemlf) = Intermediate array of reservoir for
    !                           Surface Runoffs per Gridbox (new := nmemlf = 1)
    !
    ! fgmem = Array of linear baseflow reservoir (intermediate content)
    !         At initialization it has the unit [m^3/s]
    !  friv = Array of mean riverflow = Mean Inflow per Gridbox
    !
    ! **** Other Arrays
    !
    ! area(jb) = Array of gridbox arreas, Unit = [m^2]
    !
    !
    ! **** Indices
    !
    !    jl = Longitudinal index
    !    jb = Latitudinal index
    !    il = relative change in longitude for the routing
    !    ib = relative change in latitude for the routing
    ! jlnew = jl+il
    ! jbnew = jb+ib
    !
    ! ***** Parameter and arrays for atmosphere ocean grid
    !
    !  nlon = Longitudes of atmosphere grid
    !  ngl  = Latitudes of atmosphere grid
    !
    !  oclorg = Longitudinal origin of global ocean grid
    !  ocborg = Latitudinal origin of global ocean grid
    !  ocscal = resolution  = Latitudinal width of Ocean gridbox in degree
    !
    !  aros   = atmospheric runoff array
    !  adrain = atmospheric drainage array
    !
    !  slm   = Land Sea Mask on atmosphere grid
    !  disch = Inflow array on atmospere grid
    !  xresi = Residuum (Runoff+Drainage), which results from different 
    !          land sea masks of the atmospheric grid and the 0.5 degree grid.
    !          In the latest HD model version, a start value is passed to the
    !          HD model that may include further Residual water terms that 
    !          should distributed with the discharge to close the water 
    !          balance in the coupled atmosphere ocean system.
    !

    ! mz_bk_20120730+
    ! REAL(dp) :: foslm(nlon,ngl)
    ! mz_bk_20120730-
    REAL(dp) :: zres(nlon,ngl)
    REAL(dp) :: xresi
    ! mz_bk_20120730+
    REAL(dp) :: fdata(nl, nb)
    REAL(dp) :: finp(nl, nb), finpgmem_sum, agf_k_scale
    REAL(dp) :: friv(nl,nb)
    ! REAL(dp) :: fglat(ngl)
    ! mz_bk_20120730-

    !  Parameter and switches

    INTEGER :: istep
    ! mz_bk_20120730+
    ! INTEGER :: iflow
    ! mz_bk_20120730-
    INTEGER :: jl, il, jlnew, jb, ib, jbnew
    INTEGER :: i,j,k
    ! mz_bk_20120730+
    ! , isolog, ique
    ! INTEGER :: jg, nlonp1, nglp2, isub, idum
    INTEGER :: jg, isub, idir
    ! mz_bk_20120730-

    ! mz_bk_20120730+
    ! REAL(dp) :: oclorg, ocborg, ocscal, fb, fl
    ! Ocean Grid characteristics (as used in ECHAM-Grids)
    ! Origin coordinate (usually upper left corner) & resolution
    ! Grid box centre at Longitude, Northern Border at latitude
    REAL(dp) :: ocscal
    ! mz_bk_20120730-

    REAL(dp), POINTER :: gl(:,:)

    ! mz_bk_20120730+
    ! INTEGER :: iunit
    ! LOGICAL :: lex
    ! mz_bk_20120730-

    REAL(dp) ::  zrmean

    ! NEEDED FOR WATER FLUX ADJUSTEMENT:
    ! based on mo_couple of ECHAM5 (see cosmos-1.0.0)
    REAL(dp) :: zpmebz(ngl)
    REAL(dp) :: zpmeb_glob
    REAL(dp) :: zslmz(ngl)
    REAL(dp) :: slm_glob

    ! mz_bk_20120730+
    ! IF(lhd_que) THEN
    !   ique = 1
    ! ELSE
    !   ique = 0
    ! END IF
    ! mz_bk_20120730-

    ! local scalars

  l_trig_hd         = event_state(ev_trig_hd,         current_date) 
  l_trig_hd_glacier = event_state(ev_trig_hd_glacier, current_date) 
  ! should be moved.... it is constant!
  hd_interval = event_state(ev_trig_hd)
  hd_glacier_interval = event_state(ev_trig_hd_glacier)

!  IF (ique == 1) THEN
!     write(*,*) 'HD calculation?',l_trig_hd
!      WRITE(message_text,*) 'HD calculation=',l_trig_hd
!      CALL message('hydrology_model', message_text)
!  ENDIF

  ! accumulate variables for HD-model --> done in physc
  acc_aros(:,:)    = acc_aros(:,:)+aros_2d(:,:)
  acc_adrain(:,:)  = acc_adrain(:,:)+adrain_2d(:,:)
  ! from m [apmecal_2d] to m/s [apmecal]. Done in global end
  apmecal(:,:)     = apmecal(:,:)+apmecal_2d(:,:)
  IF (lcouple) THEN
   awfre_acc(:,:)   = awfre_acc(:,:)+awfre(:,:)*time_step_len
   aifre_acc(:,:)   = aifre_acc(:,:)+aifre(:,:)*time_step_len
  ENDIF


  IF (l_trig_hd) THEN !hd_time_step 

    ! make means before transfering to HD-Model [m/s]
    zrmean = hd_interval
    !zrmean = 60._dp*60._dp*24._dp !seconds
    IF (zrmean > 0.0_dp) zrmean = 1.0_dp/zrmean
    aros(:,:)    = acc_aros(:,:)*zrmean
    adrain(:,:)  = acc_adrain(:,:)*zrmean
    IF (lcouple) THEN
      awfre_acc(:,:)   = awfre_acc(:,:)*zrmean
      aifre_acc(:,:)   = aifre_acc(:,:)*zrmean
    ENDIF
    ! mz_bk_20120928+
    IF (lnew_glac) THEN
       apmecal(:,:) = apmecal(:,:)*zrmean
    END IF
    ! mz_bk_20120928-


    ! gather data from different nodes

    gl => gl_aros
!    CALL gather_gp (gl, aros, dcg)
    CALL gather_field(gl, aros)
    gl => gl_adrain
!    CALL gather_gp (gl, adrain, dcg)
    CALL gather_field(gl, adrain)
    gl => gl_slm
!    CALL gather_gp (gl, slm, dcg)
    CALL gather_field(gl, slm)
    gl => gl_slf
!    CALL gather_gp (gl, slf, dcg)
    CALL gather_field(gl, slf)
    gl => gl_alake
!    CALL gather_gp (gl, alake, dcg)
    CALL gather_field(gl, alake)
    IF (lcouple) THEN
      ! these are accumulated!
      gl => gl_awfre
!      CALL gather_gp (gl, awfre_acc, dcg)
      CALL gather_field(gl, awfre_acc)
      gl => gl_aifre
!      CALL gather_gp (gl, aifre_acc, dcg)
      CALL gather_field(gl, aifre_acc)
    ENDIF
    ! mz_bk_20120928+
    IF (lnew_glac) THEN
       gl => gl_apmecal
       ! CALL gather_gp (gl, apmecal, dcg)
       CALL gather_field(gl, apmecal)
       gl => gl_glac
       ! CALL gather_gp (gl, glac, dcg)
       CALL gather_field(gl, glac)
    END IF
    ! mz_bk_20120928-

    ! from now on only work on IO node ... 

    IF (p_pe == p_io) THEN

       WRITE(*,*) 'HD calculations'

      ! Compute residual P-E difference over the ocean 
      ! between atmosphere and ocean model.

      xresi = 0.0_dp

      IF (lcouple) THEN
         !
         ! P-E atm. minus P-E oce.
         !
         ! mz_bk_20120730+
         ! Now done later... (mo_couple.f90)
         ! zres(:,:) = (gl_awfre(:,:)+gl_aifre(:,:))*(1.0_dp-gl_slm(:,:)) &
         !            -(gl_awfre(:,:)+gl_aifre(:,:))*(1.0_dp-gl_slf(:,:))
         ! zres(:,:) = zres(:,:)*(1.0_dp-gl_alake(:,:))                   &
         !            +(gl_awfre(:,:)+gl_aifre(:,:))*gl_alake(:,:)
         !
         !     P-E residual difference is weighted with lake mask
         !
         zres(:,:) = (gl_awfre(:,:)+gl_aifre(:,:))*gl_alake(:,:)
         ! mz_bk_20120730-
         DO jl = 1, nlon
            xresi = xresi+SUM(zres(jl,:)*gridarea(:))
         END DO

         WRITE(*,*) &
              'P-E residual difference is ', xresi, ' m**3/s = ', &
              xresi*365*86400*1.e-9_dp, ' km**3/a'
      END IF


      ! mz_bk_20120730+
      ! ! Ocean Grid characteristics (as used in ECHAM-Grids)
      ! ! Origin coordinate (usually upper left corner) & resolution
      ! ! Grid box centre at Longitude, Northern Border at latitude

      ! oclorg =   0.0_dp
      ! ocborg =  90.0_dp
      ! ocscal = 360.0_dp/nlon

      ! ! Land sea mask of Atmosphere with ocean land-sea distribution =
      ! ! Land sea mask atmosphere + lake mask

      ! foslm(:,:) = gl_slf(:,:)+gl_alake(:,:)
      ! DO jg = 1, ngl
      !   DO jl = 1, nlon
      !     IF (gl_slm(jl,jg) > 0.5_dp .AND. foslm(jl,jg) < 0.5_dp) THEN
      !       foslm(jl,jg) = gl_slm(jl,jg)
      !     END IF
      !   END DO
      ! END DO

      ocscal = fullcirc / REAL(nlon, dp)

      ! mz_bk_20120928+
      IF (lnew_glac) THEN
         !! Uwe Mikolajewicz, 2009/10/27
         !! Put P-E on glaciers into surface runoff field
         !! disable de facto the old glacier calving by setting input field
         !! to 0!!!
         !! requires ability of the HD model to transport negative runoff,
         !! which is given in ECHAM5.
         DO jg = 1, ngl
            DO jl = 1, nlon
               IF (gl_glac(jl,jg) .GT. 0.5_dp)THEN
                  gl_aros(jl,jg) = gl_aros(jl,jg) + gl_apmecal(jl,jg)
                  ! mz_bk_20120928+
                  ! has to be done on ALL processes (moved before scatter)
                  ! apmecal(jl,jg) = 0.0_dp
                  ! mz_bk_20120928-
               END IF
            END DO
         END DO
      END IF
      ! mz_bk_20120928-

      ! mz_bk_20120730-

      istep  = get_time_step()

      ! mz_bk_20120730+
      ! ! Gaussian latitudes in degrees

      ! fglat(:) = philat(:)

      ! !WRITE(message_text,*) 'fglat= ', fglat
      ! !CALL message('hydrology_model', message_text)
      ! !WRITE(message_text,*) 'gridarea= ', gridarea(1:ngl)
      ! !CALL message('hydrology_model', message_text)

      ! !  Input Runoff and simulation  of overland flow per Gridbox
      ! !  At this point in the program, it is outflow from the Gridbox

      ! iflow = 1

      ! nlonp1 = nlon+1
      ! nglp2  = ngl+2

      ! ! practically put gl_aros in the finp array (interpolated)
      ! CALL hydrology_echam(nlon, ngl, nlonp1, nglp2, gl_aros, finp, ique)
      !  Input Runoff and simulation  of overland flow per Gridbox
      !  At this point in the program, it is outflow from the Gridbox

      CALL hydrology_echam(gl_aros, finp, lhd_que)

      ! mz_bk_20120730-

      ! mz_bk_20120730+
      ! IF (iocean /= 0) THEN
      !   ! practically put gl_aros in the finp array (interpolated)
      !   ! here we correct finp on 0.5 degree grid 
      !   CALL hydrology_corr(nlon, ngl, gl_aros, gl_slm, gridarea, fglat, &
      !        finp, flag, area, xresi,                                    &
      !        oclorg, ocscal, florg, fborg, fscal, ique)
      ! ENDIF
      IF (locean) THEN
        CALL hydrology_corr(nlon, ngl, gl_aros, gl_slm, gridarea, philat, &
             finp, hd_lsm, area, xresi,                                    &
             oclorg, ocscal, florg, fborg, fscal, lhd_que)
      ENDIF
      ! mz_bk_20120730-

      !  Attention: Runoff in m/s --> Trafo with  AREA to m^3/s

      ! mz_bk_20120730+
      ! DO jb = 1, nb
      !   DO jl = 1, nl
      !     finp(jl,jb) = finp(jl,jb)*area(jb)
      !   ENDDO
      ! ENDDO
      FORALL (jb = 1:nb, jl = 1:nl)
        finp(jl,jb) = finp(jl,jb)*area(jb)
      END FORALL
      ! mz_bk_20120730-

      fdata(:,:) = 0.0_dp
 
      ! MODEL KERNEL! here the calculation are done -> output fdata,flfmem
      ! mz_bk_20120730+
      ! CALL kasglob(finp, fdata, alf_k, alf_n, iflow, flfmem, mm)
      CALL kasglob(finp, fdata, alf_k, alf_n, overlandflow, flfmem, alf_n_kas)
      ! mz_bk_20120730-
      !
      !  Reading Drainage and Computing Baseflow, Intermed. content stored in FINP

      IF (ibase /= 0) THEN

         ! mz_bk_20120730+
         ! CALL hydrology_echam(nlon, ngl, nlonp1, nglp2, gl_adrain, finp, ique)
         CALL hydrology_echam(gl_adrain, finp, lhd_que)
         ! mz_bk_20120730-

         ! mz_bk_20120730+
         ! IF (iocean /= 0) THEN
         !   CALL hydrology_corr(nlon, ngl, gl_adrain, gl_slm, gridarea, fglat, &
         !        finp, flag, area, xresi,                                    &
         !        oclorg, ocscal, florg, fborg, fscal, ique)

        ! END IF
         IF (locean) THEN
            CALL hydrology_corr(nlon, ngl, gl_adrain, gl_slm, gridarea        &
                 &            , philat                                        &
                 &            , finp, hd_lsm, area, xresi                     &
                 &            , oclorg, ocscal, florg, fborg, fscal, lhd_que)
         END IF
         ! mz_bk_20120730-

         ! *** Attention: Drainage in m/s --> Trafo with AREA to m^3/s
         ! ***    only for land points !!

         DO jl = 1, nl
            ! mz_bk_20120730+
            ! finp(jl,:) = finp(jl,:)*area(:)*flag(jl,:)
            finp(jl,:) = finp(jl,:)*area(:)*hd_lsm(jl,:)
            ! mz_bk_20120730-
         ENDDO

         ! *** Linear reservoir - Application to baseflow as done in kasglob
         ! *** the Intermediate content will be used in [m^3/s], in order to
         ! *** avoid back and forth multiplication with Unit 1 day = 86400 sec.

         ! mz_bk_20120730+
         ! fgmem(:,:) = fgmem(:,:)+finp(:,:)
         ! finp(:,:) = fgmem(:,:)/(agf_k(:,:)+1)
         ! fgmem(:,:) = fgmem(:,:)-finp(:,:)
         ! finp(:,:) = fdata(:,:)+finp(:,:)
         DO ib = 1, nb
            DO il = 1, nl
               finpgmem_sum = fgmem(il, ib) + finp(il, ib)
               agf_k_scale = finpgmem_sum / (agf_k(il, ib) + 1)
               fgmem(il, ib) = finpgmem_sum - agf_k_scale
               finp(il, ib) = fdata(il, ib) + agf_k_scale
            END DO
         END DO
         ! mz_bk_20120730-
      ELSE
         finp(:,:) = fdata(:,:)
      ENDIF

      ! ** Computing Riverflow with Input FINFL from preceeding Sub-time step

      ! mz_bk_20120730+
      ! iflow = 2
      ! mz_bk_20120730-
      friv(:,:) = 0.0_dp

      !  Computation of riverflow in MM Sub (internal) time steps
      !  i.e.  dt = 1/MM days instead of 1 day
      !
      !  Up to now a daily call of routine is forseen: may be changed to 6 hourly
      !  in later applications


      DO isub = 1, mm

         ! mz_bk_20120730+
         ! CALL kasglob(finfl, fdata, arf_k, arf_n, iflow, frfmem, mm)
         CALL kasglob(finfl, fdata, arf_k, arf_n, riverflow, frfmem, arf_n_kas)
         ! mz_bk_20120730-

        !*** Adding the riverflow
        !*** and Nullifying FINFL

        finfl(:,:) = 0.0_dp

        !*** Routing of outflow to FINFL ==> New Inflow per Gridbox

        DO jb = 1, nb
          DO jl = 1, nl
             ! mz_bk_20120730+

             ! !  *** IL, IB = relative Dircetion coordinates
             ! !  *** The 0.001-Summanden are necessary due to Rounding uncertainties

             ! ib = -(INT((fdir(jl,jb)-1)/3.0_dp+0.001_dp)-1)
             ! il = INT(((fdir(jl,jb)+2)/3.0_dp                         &
             !      -INT((fdir(jl,jb)+2)/3.0_dp+0.001_dp))*3+0.001_dp)-1

             ! !  *** Ocean point ==> FDIR = 0 ==> Pay regard by IL, IB =0

             ! idum = 1
             ! IF (fdir(jl,jb) <= 0.1_dp) idum = 0
             ! jlnew = jl+il*idum
             ! jbnew = jb+ib*idum

             ! !  *** Greenwich meridian is a boundary

             ! IF (jlnew == 0) jlnew = nl
             ! IF (jlnew == nl+1) jlnew = 1

             ! !  *** Inflow per Gridbox = Inflow+Overlandf.+Basef.+act.riverf.

             ! finfl(jlnew,jbnew) = finfl(jlnew,jbnew)+finp(jl,jb)+fdata(jl,jb)

             !  *** IL, IB = relative Dircetion coordinates
             idir = fdir(jl, jb)
             !  *** Ocean point (fdir = -1), coast(fdir = 0) or
             !      internal discharge (fdir = 5)
             ! ==> Adjust for by setting IL = IB = 0
             IF (idir > 0) THEN
                ib = 1 - (idir - 1)/3
                il = MOD(idir - 1, 3) - 1
                jlnew = MOD(jl + il - 1 + nl, nl) + 1
                jbnew = jb + ib
             ELSE
                jlnew = jl
                jbnew = jb
             END IF

             ! ib and il are guaranteed element of {-1,0,1} in this
             ! reimplementation
             !  *** Greenwich meridian is a boundary
             !  *** Inflow per Gridbox = Inflow+Overlandf.+Basef.+act.riverf.
             finfl(jlnew, jbnew) = finfl(jlnew, jbnew) &
                  + finp(jl, jb) + fdata(jl, jb)
             ! mz_bk_20120730-
          ENDDO
        ENDDO

        friv(:,:) = friv(:,:)+finfl(:,:)

        !  End of loop over Sub time steps

      ENDDO

      friv(:,:) = friv(:,:)/REAL(mm,dp)

!      IF (ique == 1) THEN
!        iunit = find_next_free_unit(90,99)
!        INQUIRE (file='discharge_0.5x0.5_diagnostics.dat', &
!             exist=lex)
!        IF (lex) THEN
!          OPEN (unit=iunit, file='discharge_0.5x0.5_diagnostics.dat', &
!               status='OLD',position='APPEND') 
!        ELSE
!          OPEN (unit=iunit, file='discharge_0.5x0.5_diagnostics.dat', &
!               status='NEW') 
!        END IF
!        WRITE(iunit) istep
!        WRITE(iunit) friv
!        CLOSE (iunit)
!      ENDIF

      !  Back trafo of Inflow to Ocean (Atmospheric) Grid

      ! mz_bk_20120730+
      ! IF (iocean /= 0) THEN
      !   CALL hydrology_to_ocean(nlon, ngl, oclorg, ocscal, fglat,  &
      !        friv, fdir, gl_disch, foslm, xresi, ique )
      ! ENDIF
      IF (locean) THEN
         CALL hydrology_to_ocean(nlon, ngl, oclorg, ocscal, philat,  &
              &                  friv, fdir, gl_disch, xresi, lhd_que)
      ENDIF
      ! mz_bk_20120730-

      ! Prepare discharge for ocean model (HOPE-C) !!!!!!!!!
      ! (should later be done outside the HD Model!)
      !
      ! Convert discharge from m**3/s to m/s for ocean model

      IF(lcouple .AND. nn == 31) THEN
          gl_disch(81,12)=gl_disch(81,12)+gl_disch(79,8)
          gl_disch(79,8)=0.0_dp
      END IF
      DO jg = 1, ngl
        DO jl = 1, nlon
          gl_disch(jl,jg) = gl_disch(jl,jg)/gridarea(jg)
        END DO
      END DO

   ENDIF !p_pe ==p_io

   ! mz_bk_20120928+
   IF (lnew_glac) apmecal(:,:) = 0.0_dp
   ! mz_bk_20120928-

   gl => gl_aros
   CALL scatter_gp (gl, aros, dcg)
   gl => gl_adrain
   CALL scatter_gp (gl, adrain, dcg)
   gl => gl_disch
   CALL scatter_gp (gl, disch_2d, dcg)
   ! mz_bk_20121001+
   ! scatter for output, in "old" version done in glacier_to_ocean,
   ! which is not called in "new" version
   ! zcalv is not used in "new" version, so set to 0.
   IF (lnew_glac) THEN
      gl => gl_disch
      CALL scatter_gp (gl, disch, dcg)
      gl_zcalv(:,:) = 0.0_dp
      gl => gl_zcalv
      CALL scatter_gp (gl, zcalv, dcg)
   END IF
   ! mz_bk_20121001-
   ! mz_bk_20120928+
!   IF (lnew_glac) THEN
!      gl => gl_apmecal
!      CALL scatter_gp (gl, apmecal, dcg)
!   END IF
   ! mz_bk_20120923-

   disch_m3s(:,:)=disch_2d(:,:)*gboxarea_2d(:,:)

   ! for safety reasons...
   ! gl_disch is later again used
   gl_disch(:,:)  = 0.0_dp
   gl_aros(:,:)   = 0.0_dp
   gl_adrain(:,:) = 0.0_dp

   ! set accumulated runoff variables zero after HD/coupling time step
   acc_aros(:,:)    = 0.0_dp
   acc_adrain(:,:)  = 0.0_dp
   IF (lcouple) THEN
     awfre_acc(:,:)   = 0.0_dp
     aifre_acc(:,:)   = 0.0_dp
   ENDIF

  ENDIF !(l_trig_hd)!hd_time_step 

  gl => gl_slm
!  CALL gather_gp (gl, slm, dcg)
  CALL gather_field(gl, slm)
  gl => gl_slf
!  CALL gather_gp (gl, slf, dcg)
  CALL gather_field(gl, slf)
  gl => gl_alake
!  CALL gather_gp (gl, alake, dcg)
  CALL gather_field(gl, alake)
  IF (lcouple) THEN
   ! these are not accumulated!!!!!!
    gl => gl_awfre
!    CALL gather_gp (gl, awfre_2d, dcg)
    CALL gather_field(gl, awfre)
    gl => gl_aifre
!    CALL gather_gp (gl, aifre_2d, dcg)
    CALL gather_field(gl, aifre)
    gl => gl_awhea
!    CALL gather_gp (gl, awhea_2d, dcg)
    CALL gather_field(gl, awhea)
    gl => gl_aicon
!    CALL gather_gp (gl, aicon_2d, dcg)
!mz_ap_20161006+ remove of field in physc.f90
    CALL gather_field(gl, ahfice)
  END IF

  IF (lcouple.and.lwater_corr) THEN
     
    ALLOCATE(qtold(kproma,ngpblks))
    ALLOCATE(qtnew(kproma,ngpblks))
    DO j=1,ngpblks
      kproma = nproma
      IF (j==ngpblks) kproma = npromz
      DO i=1,kproma
        qtold(i,j)=qtnew(i,j)
        qtnew(i,j) = 0.0_dp
        DO k=1,nlev
           qtnew(i,j)=qtnew(i,j)+(qm1(i,k,j)              &
                         +xlm1(i,k,j)+xim1(i,k,j))          &
                         *(pressi_3d(i,k+1,j)-pressi_3d(i,k,j))/g
        ENDDO
        apmebco(i,j)   =                    &
             -(qtnew(i,j)-qtold(i,j))    &
             -(rsfl_2d(i,j)+ssfl_2d(i,j)+rsfc_2d(i,j)+ssfc_2d(i,j) +qflux(i,j))*delta_time
      ENDDO
    ENDDO
    DEALLOCATE(qtold)
    DEALLOCATE(qtnew)
    gl => zpmeb
!    CALL gather_gp(gl, apmebco_2d, dcg)
    CALL gather_field(gl, apmebco)
  ENDIF


  !--------------------------------------------------------------------
  !every EVERY HD_GLACIER_TRIGGER....!!!
  !--------------------------------------------------------------------
  IF (.NOT. lnew_glac) THEN ! mz_bk_20120928
     IF (l_trig_hd_glacier) THEN !hd_glacier_time_step 
        ! to be averaged
        apmecal(:,:) = apmecal(:,:)/hd_glacier_interval
        gl => gl_apmecal
        !    CALL gather_gp (gl, apmecal, dcg) 
        CALL gather_field(gl, apmecal)

        gl => gl_glac
        !    CALL gather_gp (gl, glac, dcg)
        CALL gather_field(gl, glac)

        gl => gl_disch
        !    CALL gather_gp (gl, disch_2d, dcg)
        CALL gather_field(gl, disch_2d)

        IF (p_pe ==p_io) THEN  
           ! mz_bk_20120730+
           ! IF (ique == 1) THEN
           IF (lhd_que) THEN
              ! mz_bk_20120730-
              write(*,*) 'glacier to ocean (calvin model)'
           ENDIF
           ! here we calculate gl_zcalv....
           ! mz_bk_20120730+
           ! CALL glacier_to_ocean
           CALL glacier_to_ocean(lhd_que)
           ! mz_bk_20120730-
           DO jg = 1, ngl
              DO jl = 1, nlon
                 gl_disch(jl,jg) = gl_disch(jl,jg)+gl_zcalv(jl,jg)/gridarea(jg)
              ENDDO
           ENDDO
        ENDIF !p_pe ==p_io

        apmecal(:,:) = 0.0_dp

        gl => gl_disch
        CALL scatter_gp (gl, disch, dcg)
        gl => gl_zcalv
        CALL scatter_gp (gl, zcalv, dcg)

     ENDIF !hd_glacier_time_step
  END IF ! mz_bk_20120928 ! .not. lnew_glac

  !--------------------------------------------------------------------
  !every TIME STEP!!!
  !--------------------------------------------------------------------
  gl => gl_disch
!  CALL gather_gp (gl, disch_2d, dcg)
  CALL gather_field(gl, disch_2d)
  ! mz_bk_20120423+
  gl => gl_zcalv
!  CALL gather_gp (gl, zcalv, dcg)
  CALL gather_field(gl, zcalv)
  ! mz_bk_20120423-

  IF (p_pe ==p_io) THEN  
     ! mz_bk_20120730+
     ! IF (ique == 1) THEN
     IF (lhd_que) THEN
     ! mz_bk_20120730-
      write(*,*) 'Add discharge to fresh water flux for ocean'
    ENDIF
    ! Add discharge to fresh water flux for ocean
    IF (lcouple) THEN
       DO jg = 1, ngl
          DO jl = 1, nlon
             IF (gl_slf(jl,jg)+gl_alake(jl,jg) < 0.95_dp) THEN
                gl_awfre(jl,jg) = gl_awfre(jl,jg)+gl_disch(jl,jg) &
                     /(1.0_dp-gl_slf(jl,jg)-gl_alake(jl,jg))
             END IF
             IF (.NOT. lnew_glac) THEN ! mz_bk_20120928
                ! Add calved water to fresh water flux into the ocean
                ! and subtract latent heats from heat flux and conductive heat
                gl_awfre(jl,jg) = gl_awfre(jl,jg)+gl_zcalv(jl,jg)/gridarea(jg) &
                     /MAX(1.0_dp-gl_slf(jl,jg)-gl_alake(jl,jg),1.e-6_dp)
                ! mz_bk_20120730+
                ! gl_awhea(jl,jg) = gl_awhea(jl,jg)-gl_zcalv(jl,jg)/gridarea(jg) &
                !      *rhoh2o*alf  
                ! gl_aicon(jl,jg) = gl_aicon(jl,jg)-gl_zcalv(jl,jg)/gridarea(jg) &
                !      *rhoh2o*alf
                gl_awhea(jl,jg) = gl_awhea(jl,jg)-gl_zcalv(jl,jg)/gridarea(jg) &
                     /MAX(1.0_dp-gl_slf(jl,jg)-gl_alake(jl,jg),1.e-6_dp)*rhoh2o*alf
                gl_aicon(jl,jg) = gl_aicon(jl,jg)-gl_zcalv(jl,jg)/gridarea(jg) &
                     /MAX(1.0_dp-gl_slf(jl,jg)-gl_alake(jl,jg),1.e-6_dp)*rhoh2o*alf
                ! mz_bk_20120730-
             END IF ! mz_bk_20120928
          END DO
       END DO
    ENDIF

    ! HERE WE ADD ALSO THE WATER CORRECTION OF COUPLED MODEL
    ! done every time step and not at the coupling time!!!
    ! we put it here because most of the variables needed
    ! were already calculated...
    IF (lwater_corr.and.lcouple) THEN

      ! as done in ioinitial.f90
      DO jg=1,ngl
        zslmz(jg)=SUM(gl_slm(1:nlon,jg)                            &
                     +gl_alake(1:nlon,jg))*gl_budw(jg)
      END DO
      slm_glob=SUM(zslmz)
      !CALL p_bcast (slm_glob,  p_io)

      ! global mean
      DO jg=1,ngl
        zpmebz(jg)=SUM(zpmeb(1:nlon,jg))*gl_budw(jg)
      END DO

      zpmeb_glob=SUM(zpmebz)

      ! Apply global p-e correction to fresh water flux

      zcorr=zpmeb_glob/((1.-slm_glob)*rho_H2O*delta_time)

      ! not more needed to mask 
      ! (done in messy_a2o_e5.f90, during transformation)
      WHERE (gl_slf(:,:).lt.0.5_dp.and.gl_alake(:,:).lt.0.5_dp)
      gl_awfre(:,:) = gl_awfre(:,:) + zcorr
      END WHERE

    ENDIF

  ENDIF !p_pe ==p_io
  !----------------------------------------------------------------

  IF (lcouple) THEN
    gl => gl_awfre
    CALL scatter_gp (gl, awfre, dcg)
    gl => gl_awhea
    CALL scatter_gp (gl, awhea, dcg)
    gl => gl_aicon
    CALL scatter_gp (gl, ahfice, dcg)
  END IF


  END SUBROUTINE hd_global_end

  ! ---------------------------------------------------------------------------

  SUBROUTINE hd_free_memory

    ! mz_bk_20120730+
    ! CALL cleanup_hydrology_slm_invariants
    CALL cleanup_hydrology_slm_invar
    ! mz_bk_20120730-

    IF (p_pe == p_io) THEN

      DEALLOCATE  (alf_k)
      DEALLOCATE  (alf_n)
      DEALLOCATE  (arf_k)
      DEALLOCATE  (arf_n)
      DEALLOCATE  (agf_k)
      DEALLOCATE  (fdir)
      ! mz_bk_20120730+
      ! DEALLOCATE  (flag)
      DEALLOCATE  (hd_lsm)
      ! mz_bk_20120730-
      ! mz_bk_20120730+
      ! DEALLOCATE  (friv)
      ! DEALLOCATE  (finp)
      ! DEALLOCATE  (fdata)
      ! mz_bk_20120730-
      DEALLOCATE  (area)
      ! mz_bk_20120730+
      DEALLOCATE  (oclook_cache)
      DEALLOCATE  (intpol_mapping)
      ! mz_bk_20120730-
      
    END IF

    DEALLOCATE (gl_aros)
    DEALLOCATE (gl_adrain)     
    DEALLOCATE (gl_disch)  
    DEALLOCATE (gl_zcalv)  
    DEALLOCATE (gl_slm)    
    DEALLOCATE (gl_slf)    
    DEALLOCATE (gl_alake)  
    IF (lcouple) THEN
      DEALLOCATE (gl_awfre)  
      DEALLOCATE (gl_aifre)  
      DEALLOCATE (gl_awhea)  
    ENDIF
    DEALLOCATE (gl_aicon)  
    DEALLOCATE (gl_apmecal)  
    DEALLOCATE (gl_glac)  
!mz_ap_20090306+
    DEALLOCATE (zpmeb)

  END SUBROUTINE hd_free_memory

  ! ---------------------------------------------------------------------------

! ===========================================================================
! PRIVATE hd INTERFACE ROUTINES
! ===========================================================================

  SUBROUTINE hd_read_start

    USE messy_main_timer,             ONLY: lstart
    USE mo_io

    TYPE (FILE_INFO)  :: fileinfo
    ! files reading 
    LOGICAL :: lex
    ! File names
    INTEGER nvarid, fileid, i
    CHARACTER(len= 7) :: varname
    CHARACTER(len=100) :: yodnpar, yodnres
    ! mz_bk_20120730+
    REAL(dp) :: hd_dp_read(nl,nb)
    ! mz_bk_20120730-

    yodnpar = TRIM(inithd_files_path)//'/hdpara.nc'
    
    ! Read parameter: Land sea mask, RDF, ...
    
    INQUIRE (file=yodnpar, exist=lex)
    IF (.NOT. lex) THEN
      WRITE (*,*) 'Could not open file <',TRIM(yodnpar),'>'
      CALL finish ('read_hydrology', 'run terminated.')
    ENDIF
    
    ! reading grid file....

    fileinfo%opened = .FALSE.
    CALL IO_open (yodnpar, fileinfo, IO_READ)
    CALL message('read_hydrology',  'Reading hdpara from file '//TRIM(yodnpar))
    
    fileID = fileinfo%file_id
    
    CALL IO_inq_varid (fileID, 'FLAG', nvarid)
    ! mz_bk_20120730+
    ! CALL IO_get_var_double (fileID, nvarid, flag)
    CALL IO_get_var_double (fileID, nvarid, hd_dp_read)
    hd_lsm = REAL(hd_dp_read, sp)
    ! mz_bk_20120730-
    CALL IO_inq_varid (fileID, 'FDIR', nvarid)
    ! mz_bk_20120730+
    ! CALL IO_get_var_double (fileID, nvarid, fdir)
    CALL IO_get_var_double (fileID, nvarid, hd_dp_read)
    ! convert back to integer directions pointlessly stored as double
    fdir = MIN(INT(hd_dp_read), 9)
    ! mz_bk_20120730-
    CALL IO_inq_varid (fileID, 'ALF_K', nvarid)
    CALL IO_get_var_double (fileID, nvarid, alf_k)
    CALL IO_inq_varid (fileID, 'ALF_N', nvarid)
    CALL IO_get_var_double (fileID, nvarid, alf_n)
    CALL IO_inq_varid (fileID, 'ARF_K', nvarid)
    CALL IO_get_var_double (fileID, nvarid, arf_k)
    CALL IO_inq_varid (fileID, 'ARF_N', nvarid)
    CALL IO_get_var_double (fileID, nvarid, arf_n)
    CALL IO_inq_varid (fileID, 'AGF_K', nvarid)
    CALL IO_get_var_double (fileID, nvarid, agf_k)
    
    CALL IO_close(fileinfo)

    ! reading starting conditions (not needed at restart!)

    yodnres = TRIM(inithd_files_path)//'/hdstart.nc'

    if (lstart) THEN
    ! in case of rerun we use the rerun information!
      INQUIRE (file=yodnres, exist=lex)
      IF (.NOT. lex) THEN
        WRITE (*,*) 'Could not open file <',TRIM(yodnpar),'>'
        CALL finish ('read_hydrology', 'run terminated.')
      ENDIF
    
      ! reading initial conditions....

      fileinfo%opened = .FALSE.
      CALL IO_open (yodnres, fileinfo, IO_READ)
      CALL message('read_hydrology',  'Reading hdstart from file '//TRIM(yodnres))
    
      fileID = fileinfo%file_id

      CALL IO_inq_varid (fileID, 'FLFMEM', nvarid)
      CALL IO_get_var_double (fileID, nvarid, flfmem)
      
      varname = 'FRFMEM'
      DO i=1, nmemrf
        WRITE(varname(7:7), '(i1)') i
        CALL IO_inq_varid (fileID, varname, nvarid)
        CALL IO_get_var_double (fileID, nvarid, frfmem(:,:,I))
      ENDDO
      
      CALL IO_inq_varid (fileID, 'FGMEM', nvarid)
      CALL IO_get_var_double (fileID, nvarid, fgmem)
      
      ! mz_bk_20120730+
      ! fgmem(:,:) = fgmem(:,:)*flag(:,:)
      fgmem(:,:) = fgmem(:,:)*REAL(hd_lsm(:,:), dp)
      ! mz_bk_20120730-
      
      CALL IO_inq_varid (fileID, 'FINFL', nvarid)
      CALL IO_get_var_double (fileID, nvarid, finfl)
      
      CALL IO_close(fileinfo)
    
    endif !lstart

  END SUBROUTINE hd_read_start

!-------------------------------------------------------------------------------
  ! mz_bk_20120730+
  ! Land sea mask of Atmosphere with ocean land-sea distribution =
  ! Land sea mask atmosphere + lake mask
  SUBROUTINE slm_lake_join(nlat, nlon, slm_dest, slm_src, slf_src, alake_src)
    INTEGER, INTENT(in) :: nlat, nlon
    REAL(dp), INTENT(out) :: slm_dest(nlon, nlat)
    REAL(dp), INTENT(in) :: slf_src(nlon, nlat), slm_src(nlon, nlat), &
         alake_src(nlon, nlat)

    slm_dest = MERGE(slm_src, slf_src + alake_src, &
         slm_src > 0.5_dp .AND. slf_src + alake_src < 0.5_dp)
  END SUBROUTINE slm_lake_join
  ! mz_bk_20120730-

  ! mz_bk_20120730+
!   SUBROUTINE hydrology_echam (nlon, nlat, nlonp1, nlatp2, code, tocode, ique)

!     !*************************************************************************
!     !
!     ! **** This program interpolates data from Gaussian grids to a half
!     !    degree grid
!     !
!     !  Programmierung und Entwicklung: Uwe Schulzweida (echamto30min)
!     !  Modified to Subroutine by Stefan Hagemann -- September 1995
!     !
!     ! ***** Version 1.1 -- Dezember 1995
!     !          Instead of longitude centred coordinate array trcode, now
!     !          the 0.5 degree coordinate field tocode is given back to the calling
!     !          routine which has a perfect boundary with the Northpole/dateline
!     !          --> Origin has centre coordinate 89.75 N, -179.75 W
!     !
!     ! ***** Version 2.0 -- November 1999
!     !   Since the input array to be transformed is only passed to ECHREAD
!     !   but not read in ECHREAD itself, in ECHREAD only 
!     !   the transformation of the input array to 0.5 degree is done.
!     !   Therefor, new calling parameter and routine names are set/given.
!     !          alt: SUBROUTINE echread(luinp, ihead, tocode, istep, ique)
!     !          neu: SUBROUTINE hydrology_model(code, tocode, ique)
!     !
!     ! ****** List of variables
!     !
!     !  trcode = Interpolated Array
!     !  ique   = Log-output switch  ( 0 = No Log-Output )

!     USE messy_main_data_bi, ONLY: philat, philon

!     INTEGER, INTENT(in) :: nlat, nlon, nlonp1, nlatp2, ique

!     REAL(dp), PARAMETER :: umfang = 360.0_dp

!     REAL(dp) :: code(nlon,nlat), xlon(nlon), xlat(nlat)
!     REAL(dp) :: acode(nlonp1,nlatp2), axlon(nlonp1), axlat(nlatp2)
!     REAL(dp) :: trcode(nl,nb), xr(nl), yr(nb)
!     REAL(dp) :: tocode(nl,nb)

!     INTEGER :: j, jlat, jlon

!     ! definition of the input data grid

!     DO j = 1, nlat
!       xlat(j)    = philat(j)
!       axlat(j+1) = xlat(j)
!     ENDDO

!     axlat(1)      =  90.0_dp
!     axlat(nlatp2) = -90.0_dp

!     IF (ique /= 0) THEN
!       WRITE(message_text,*) xlat(1), xlat(2), xlat(nlat)
!       CALL message('hydrology_echam', message_text)
!       WRITE(message_text,*) axlat(1), axlat(2), axlat(nlatp2)
!       CALL message('hydrology_echam', message_text)
!     ENDIF

!     DO j = 1, nlon
!       xlon(j)  = philon(j)
!       axlon(j) = xlon(j)
!     ENDDO

!     axlon(nlonp1) = 360.0_dp

!     IF (ique /= 0) THEN
!       WRITE(message_text,*) xlon(1), xlon(2), xlon(nlon)
!       CALL message('hydrology_echam', message_text)
!       WRITE(message_text,*) axlon(1), axlon(2), axlon(nlonp1)
!       CALL message('hydrology_echam', message_text)
!     ENDIF

!     ! definition of the output data grid

!     DO j = 1, nl
!       xr(j) = 0.25_dp+(j-1)/REAL(nl,dp)*umfang
!     END DO

!     IF (ique /= 0) THEN
!       WRITE(message_text,*) xr(1), xr(2), xr(nl)
!       CALL message('hydrology_echam', message_text)
!     END IF

!     DO j = 1, nb
!       yr(j) = 0.25_dp-0.25_dp*umfang+0.5_dp*umfang*(j-1)/REAL(nb,dp)
!       yr(j) = -yr(j)
!     ENDDO

!     IF (ique /= 0) THEN
!       WRITE(message_text,*) yr(1), yr(2), yr(nb)
!       CALL message('hydrology_echam', message_text)
!     END IF

!     ! compare expected timestep with next one to read

!     DO jlat = 1, nlat
!       DO jlon = 1, nlon
!         acode(jlon,jlat+1) = code(jlon,jlat)
!       ENDDO
!     ENDDO

!     DO jlat = 2, nlatp2-1
!       acode(nlonp1,jlat) = acode(1,jlat)
!     ENDDO

!     DO jlon = 1, nlonp1
!       acode(jlon,1)      = acode(jlon,2)
!       acode(jlon,nlatp2) = acode(jlon,nlatp2-1)
!     ENDDO

!     ! interpolation to output grid

!     CALL intpol(nlonp1, nlatp2, acode, axlon, axlat, nl, nb, trcode, xr, yr)

!     ! transformation on bounded coordinates

!     DO jlat = 1, nb
!       DO jlon = 1, nl/2
!         tocode(jlon,jlat)      = trcode(jlon+nl/2,jlat)
!         tocode(jlon+nl/2,jlat) = trcode(jlon,jlat)
!       ENDDO
!     ENDDO

!   END SUBROUTINE hydrology_echam

! !-------------------------------------------------------------------------------

!   SUBROUTINE intpol (nxm, nym, fieldm, xm, ym, nx, ny, field, x, y)

!     INTEGER, INTENT(in) :: nxm, nym, nx, ny

!     REAL(dp), INTENT(in)  :: x(nx),y(ny)
!     REAL(dp), INTENT(in)  :: xm(nxm),ym(nym)
!     REAL(dp), INTENT(in)  :: fieldm(nxm,nym)
!     REAL(dp), INTENT(out) :: field(nx,ny)

!     INTEGER :: irun, jj, j, ii, i

!     irun = 0
!     DO jj = 2, nym
!       DO j = 1, ny
!         ! if (y(j) < ym(jj-1) .or. y(j) > ym(jj)) cycle
!         IF (y(j) < MIN(ym(jj-1),ym(jj)) .OR.   &
!             y(j) > MAX(ym(jj-1),ym(jj))) CYCLE
!         irun = irun+1
!         DO ii = 2, nxm
!           DO i = 1, nx
!             IF(x(i) < xm(ii-1) .OR. x(i) > xm(ii)) CYCLE
!             field(i,j) = fieldm(ii-1,jj-1)*(x(i)-xm(ii))*(y(j)-ym(jj)) &
!                         /((xm(ii-1)-xm(ii))*(ym(jj-1)-ym(jj)))         &
!                         +fieldm(ii,jj-1)*(x(i)-xm(ii-1))*(y(j)-ym(jj)) &
!                         /((xm(ii)-xm(ii-1))*(ym(jj-1)-ym(jj)))         &
!                         +fieldm(ii-1,jj)*(x(i)-xm(ii))*(y(j)-ym(jj-1)) &
!                         /((xm(ii-1)-xm(ii))*(ym(jj)-ym(jj-1)))         &
!                         +fieldm(ii,jj)*(x(i)-xm(ii-1))*(y(j)-ym(jj-1)) &
!                         /((xm(ii)-xm(ii-1))*(ym(jj)-ym(jj-1)))
!           ENDDO
!         ENDDO
!       ENDDO
!     ENDDO

!   END SUBROUTINE intpol

  SUBROUTINE hydrology_echam (field_in, field_out, lhd_que)

    !*************************************************************************
    !
    ! **** This program interpolates data from Gaussian grids to a half
    !    degree grid
    !
    !  Programmierung und Entwicklung: Uwe Schulzweida (echamto30min)
    !  Modified to Subroutine by Stefan Hagemann -- September 1995
    !
    ! ***** Version 1.1 -- Dezember 1995
    !          Instead of longitude centred coordinate array trcode, now
    !          the 0.5 degree coordinate field field_out is given back to the calling
    !          routine which has a perfect boundary with the Northpole/dateline
    !          --> Origin has centre coordinate 89.75 N, -179.75 W
    !
    ! ***** Version 2.0 -- November 1999
    !   Since the input array to be transformed is only passed to ECHREAD
    !   but not read in ECHREAD itself, in ECHREAD only
    !   the transformation of the input array to 0.5 degree is done.
    !   Therefor, new calling parameter and routine names are set/given.
    !          alt: SUBROUTINE echread(luinp, ihead, field_out, istep, lhd_que)
    !          neu: SUBROUTINE hydrology_model(field_in, field_out, lhd_que)
    !
    ! ****** List of variables
    !
    !  field_out = Interpolated, transposed Array
    !  lhd_que   = Log-output switch  ( 0 = No Log-Output )
    USE messy_main_grid_def_mem_bi,   ONLY: nlon, ngl 
    USE messy_main_grid_def_bi,       ONLY: philat, philon

    REAL(dp), INTENT(in) :: field_in(nlon,ngl)
    REAL(dp), INTENT(out) :: field_out(nl,nb)
    LOGICAL, INTENT(in) :: lhd_que

    REAL(dp) :: acode(nlon + 1,ngl + 2), axlon(nlon + 1), axlat(ngl + 2)
    REAL(dp) :: xr(nl), yr(nb)

    INTEGER :: jlat, jlon


    CALL intpol_coord_axis_setup(axlon, axlat, xr, yr)

    IF (lhd_que) THEN
      WRITE(message_text,*) philat(1), philat(2), philat(ngl)
      CALL message('hydrology_echam', message_text)
      WRITE(message_text,*) axlat(1), axlat(2), axlat(ngl + 2)
      CALL message('hydrology_echam', message_text)

      WRITE(message_text,*) philon(1), philon(2), philon(nlon)
      CALL message('hydrology_echam', message_text)
      WRITE(message_text,*) axlon(1), axlon(2), axlon(nlon + 1)
      CALL message('hydrology_echam', message_text)

      WRITE(message_text,*) xr(1), xr(2), xr(nl)
      CALL message('hydrology_echam', message_text)

      WRITE(message_text,*) yr(1), yr(2), yr(nb)
      CALL message('hydrology_echam', message_text)
    END IF

    ! compare expected timestep with next one to read

    DO jlat = 1, ngl
      DO jlon = 1, nlon
        acode(jlon,jlat+1) = field_in(jlon,jlat)
      ENDDO
    ENDDO

    DO jlat = 2, ngl + 1
      acode(nlon + 1,jlat) = acode(1,jlat)
    ENDDO

    DO jlon = 1, nlon + 1
      acode(jlon,1)      = acode(jlon,2)
      acode(jlon,ngl + 2) = acode(jlon,ngl + 2-1)
    ENDDO

    ! interpolation to output grid
    ! (includes transformation on bounded coordinates)
    CALL intpol_with_mapping(nlon + 1, ngl + 2, acode, axlon, axlat, &
         nl, nb, field_out, xr, yr, intpol_mapping)
  END SUBROUTINE hydrology_echam

  SUBROUTINE intpol_coord_axis_setup(axlon, axlat, xr, yr)
    USE messy_main_grid_def_mem_bi,   ONLY: nlon, ngl 
    USE messy_main_grid_def_bi,       ONLY: philat, philon
    REAL(dp), INTENT(inout) :: axlon(nlon + 1), axlat(ngl + 2), xr(nl), yr(nb)

    INTEGER :: j

    ! definition of the echam gp data grid
    axlat(1)       =  90.0_dp
    axlat(2:ngl+1) =  philat(1:ngl)
    axlat(ngl + 2) = -90.0_dp

    axlon(1:nlon)  = philon(:)
    axlon(nlon + 1) = fullcirc

    ! definition of hd data grid
    DO j = 1, nl
      xr(j) = 0.5_dp*hd_scal_lon + (j-1)/REAL(nl,dp)*fullcirc
    END DO

    DO j = 1, nb
      ! perhaps use this formulation? (tj, 20091009)
      ! yr(j) = -(0.5*hd_scal_lat + fullcirc * hd_scal_lat * ( (j-1) &
      !          / REAL(nb,dp) - 0.5))
      yr(j) = -(0.5_dp*hd_scal_lat - 0.5_dp*hd_scal_lat * fullcirc &
           &    + hd_scal_lat * fullcirc * REAL(j-1, dp) / REAL(nb,dp))
    ENDDO
  END SUBROUTINE intpol_coord_axis_setup


  !> interpolate field src to dest
  SUBROUTINE intpol (src_i_size, src_j_size, src, src_x, src_y, &
       dest_i_size, dest_j_size, dest, dest_x, dest_y)

    INTEGER, INTENT(in) :: src_i_size, src_j_size, dest_i_size, dest_j_size

    REAL(dp), INTENT(in)  :: dest_x(dest_i_size), dest_y(dest_j_size)
    REAL(dp), INTENT(in)  :: src_x(src_i_size), src_y(src_j_size)
    REAL(dp), INTENT(in)  :: src(src_i_size,src_j_size)
    REAL(dp), INTENT(out) :: dest(dest_i_size,dest_j_size)

    INTEGER :: irun, js, jd, is, id, idt

    irun = 0
    DO js = 2, src_j_size
      DO jd = 1, dest_j_size
        ! if (dest_y(jd) < src_y(js-1) .or. dest_y(jd) > src_y(js)) cycle
        IF (dest_y(jd) < MIN(src_y(js - 1), src_y(js)) .OR.   &
             dest_y(jd) > MAX(src_y(js - 1), src_y(js))) CYCLE
        irun = irun+1
        DO is = 2, src_i_size
          DO id = 1, dest_i_size
            IF(dest_x(id) < src_x(is-1) .OR. dest_x(id) > src_x(is)) CYCLE
            idt = MOD(id - 1 + dest_i_size/2, dest_i_size) + 1
            dest(idt,jd) = src(is - 1, js - 1) &
                 * (dest_x(id) - src_x(is)) * (dest_y(jd) - src_y(js)) &
                 / ((src_x(is - 1) - src_x(is)) * (src_y(js - 1) - src_y(js))) &
                 + src(is, js - 1) &
                 * (dest_x(id) - src_x(is - 1)) * (dest_y(jd) - src_y(js)) &
                 / ((src_x(is) - src_x(is - 1)) * (src_y(js - 1) - src_y(js))) &
                 + src(is - 1, js) &
                 * (dest_x(id) - src_x(is)) * (dest_y(jd) - src_y(js - 1)) &
                 / ((src_x(is - 1) - src_x(is)) * (src_y(js) - src_y(js - 1))) &
                 + src(is, js) &
                 * (dest_x(id) - src_x(is - 1)) * (dest_y(jd) - src_y(js - 1)) &
                 / ((src_x(is) - src_x(is - 1)) * (src_y(js) - src_y(js - 1)))
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE intpol

  !> find midpoint is,js for every index pair from field
  !> src(src_i_size, src_j_size) to dest(dest_i_size, dest_j_size)
  !> every element cart_idx_2d(is,js) of mapping(id,jd) later defines
  !> that dest(id,jd) will be computed from
  !> src(is,js), src(is - 1, js), src(is, js - 1), src(is - 1, js - 1)
  SUBROUTINE intpol_compute_mapping(mapping, &
       src_i_size, src_j_size, src_x, src_y, &
       dest_i_size, dest_j_size, dest_x, dest_y)

    INTEGER, INTENT(in) :: src_i_size, src_j_size, dest_i_size, dest_j_size

    REAL(dp), INTENT(in)  :: dest_x(dest_i_size),dest_y(dest_j_size)
    REAL(dp), INTENT(in)  :: src_x(src_i_size),src_y(src_j_size)
    TYPE(cart_idx_2d), INTENT(out) :: mapping(dest_i_size,dest_j_size)

    INTEGER :: js, jd, is, id

    mapping = cart_idx_2d(-1, -1)
    DO js = 2, src_j_size
      DO jd = 1, dest_j_size
        ! if (dest_y(jd) < src_y(js-1) .or. dest_y(jd) > src_y(js)) cycle
        IF (dest_y(jd) < MIN(src_y(js - 1), src_y(js)) .OR.   &
            dest_y(jd) > MAX(src_y(js - 1), src_y(js))) CYCLE
        DO is = 2, src_i_size
          DO id = 1, dest_i_size
            IF(dest_x(id) < src_x(is-1) .OR. dest_x(id) > src_x(is)) CYCLE
            mapping(id, jd) = cart_idx_2d(is, js)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE intpol_compute_mapping

  !> interpolate field src to dest
  SUBROUTINE intpol_with_mapping(src_i_size, src_j_size, src, src_x, src_y, &
       dest_i_size, dest_j_size, dest, dest_x, dest_y, mapping)

    INTEGER, INTENT(in) :: src_i_size, src_j_size, dest_i_size, dest_j_size

    REAL(dp), INTENT(in)  :: dest_x(dest_i_size),dest_y(dest_j_size)
    REAL(dp), INTENT(in)  :: src_x(src_i_size),src_y(src_j_size)
    REAL(dp), INTENT(in)  :: src(src_i_size,src_j_size)
    REAL(dp), INTENT(out) :: dest(dest_i_size,dest_j_size)
    TYPE(cart_idx_2d), INTENT(in) :: mapping(dest_i_size, dest_j_size)
    INTEGER :: js, jd, is, id, idt

    DO jd = 1, dest_j_size
      DO id = 1, dest_i_size
        IF (mapping(id,jd)%ilat /= -1) THEN
          is = mapping(id,jd)%ilon
          js = mapping(id,jd)%ilat
          idt = id + dest_i_size/2 - MERGE(dest_i_size, 0, id > dest_i_size/2)
          dest(idt,jd) = src(is - 1, js - 1) &
               * (dest_x(id) - src_x(is)) * (dest_y(jd) - src_y(js)) &
               / ((src_x(is - 1) - src_x(is)) * (src_y(js - 1) - src_y(js))) &
               + src(is, js - 1) &
               * (dest_x(id) - src_x(is - 1)) * (dest_y(jd) - src_y(js)) &
               / ((src_x(is) - src_x(is - 1)) * (src_y(js - 1) - src_y(js))) &
               + src(is - 1, js) &
               * (dest_x(id) - src_x(is)) * (dest_y(jd) - src_y(js - 1)) &
               / ((src_x(is - 1) - src_x(is)) * (src_y(js) - src_y(js - 1))) &
               + src(is, js) &
               * (dest_x(id) - src_x(is - 1)) * (dest_y(jd) - src_y(js - 1)) &
               / ((src_x(is) - src_x(is - 1)) * (src_y(js) - src_y(js - 1)))
        END IF
      END DO
    END DO
  END SUBROUTINE intpol_with_mapping

  ! mz_bk_20120730-


!-------------------------------------------------------------------------------
  ! mz_bk_20120730+
!   SUBROUTINE hydrology_corr(nlon, nlat, fatmos, foclsm, aoarea, fglat, &
!                             fdat, flag, area, xresi,                   &
!                             oclorg, ocscal, florg, fborg, fscal, ique)


!     ! **** Correction of the atmospheric grid -> 0.5 degree transformation
!     !
!     ! ***** Version 1.0 - November 1999
!     ! Programmed and developed by Stefan Hagemann, MPI
!     ! Remark: Input data FATMOS (Runoff or Drainage) should have the unit m/s.
!     !
!     ! ***** Version 1.1 - January 2001
!     ! Longitude-Index-Correction
!     !
!     ! ***** Version 2.1 - January 2001
!     ! ECHAM5- Version incl. Gaussian latitudes
!     !
!     !
!     ! **** Global Arrays:
!     !
!     !  fdat = Data array at 0.5 degree
!     !  flag = Land mask array at 0.5 degree
!     ! area(jb) = Array of Gridbox areas, Unit = [m^2]
!     !
!     ! ***** Parameter and arrays for atmosphere ocean grid
!     !
!     !    nlon = Longitudes of global ocean grid
!     !    nlat = Latitudes of global ocean grid
!     !
!     !  oclorg = Longitudinal origin of global ocean grid
!     !  ocscal = resolution/scale = Width of an ocean Gridbox in degree
!     !
!     !  fatmos = atmospheric data array
!     !  foclsm = Land Sea Mask on Ocean grid
!     !  aoarea = Area per gridbox per latitude of atmospheric grid [m^2]
!     !   fglat = Gaussian latitude of global ocean grid (centre coordinates)
!     !
!     !   xresi = Residuum (Runoff+Drainage), which results from different 
!     !           land sea masks of the atmospheric grid and the 0.5 degree grid.
!     !           In the latest HD model version, a start value is passed to the
!     !           HD model that may include further residual water terms that should
!     !           distributed with the discharge to close the water balance in the
!     !           coupled atmosphere ocean system.
!     !
!     !     ndd = If no land point is found as direct neighbour,
!     !           it is searched in NWSE direction until the maximum distance of
!     !           NDD Boxes is reached.
!     !    ique = Log-Output switch ( 0 = No Log-output to STDOUT)
!     !

!     INTEGER, PARAMETER :: ndd = 3

!     INTEGER, INTENT(in) :: nlon, nlat, ique

!     ! Input fields AO-Grid

!     REAL(dp) :: fatmos(nlon,nlat), aoarea(nlat), foclsm(nlon,nlat)
!     REAL(dp) :: fglat(nlat)

!     REAL(dp) :: fdat(nl,nb), flag(nl,nb), area(nb)
!     REAL(dp) :: x1, x2, xresi
!     REAL(dp) :: oclorg, ocscal
!     REAL(dp) :: florg, fborg, fscal

!     INTEGER :: jb, jl, jlon, jlat, idd

!     REAL(dp) :: xanf(nlat), xend(nlat)
!     REAL(dp) :: xjlon(nl), fb(nl), fl(nl), xx1(nl), xx2(nl)
!     INTEGER :: jjlat(nl), jjlon(nl)
!     LOGICAL :: lset(nl), laction(nl)


!     ! Precompute xanf and xend

!     DO jlat = 1, nlat
!       IF (jlat == 1) THEN
!         xanf(jlat) = 90.0_dp
!       ELSE
!         xanf(jlat) = (fglat(jlat-1)+fglat(jlat))*0.5_dp
!       ENDIF
!       IF (jlat == nlat) THEN
!         xend(jlat) = -90.0_dp
!       ELSE
!         xend(jlat) = (fglat(jlat+1)+fglat(jlat))*0.5_dp
!       ENDIF
!     END DO

!     ! OA = Ocean-Atmosphere Grid, HD = 0.5 Grad HD-Model Grid
    
!     DO jb = 1, nb-6
!       ! First, compute fb and fl for all jl values
!       DO jl = 1, nl
!         ! grid box centre:
!         fb(jl) = fborg-REAL(jb,dp)*fscal+0.5_dp*fscal
!         fl(jl) = REAL(jl,dp)*fscal+florg-0.5_dp*fscal
!         IF (fl(jl) >= 180.0_dp) fl(jl) = fl(jl)-360.0_dp
!       END DO
      
!       ! Corresponding Index in Ocean Grid
!       ! Latitude

! !!! without Gauss:    XJLAT = (OCBORG-FB) / OCSCAL+1.  , OCBORG = Borderline
! !!!                 JLAT = INT(XJLAT+0.00001)

!       ! and compute jlat and jlog for all jl values

!       lset = .TRUE.
!       DO jlat = 1, nlat
!         DO jl = 1, nl
!           IF(fb(jl) <= xanf(jlat) .AND. fb(jl) > xend(jlat) .AND. lset(jl)) THEN
!             jjlat(jl) = jlat
!             lset(jl)  = .FALSE.
!           ENDIF
!         ENDDO
!       END DO

!       lset=(jjlat == nlat+1)

!       IF(ANY(lset))  THEN
!          DO jl = 1, nl
!            jlat = jjlat(jl)
!            IF (jlat == nlat+1) THEN
!              WRITE(message_text,*) ' error in jlat=', jlat
!              CALL message ('hydrology_corr', message_text)
!              jjlat(jl) = nlat
!            ENDIF
!          END DO
!       END IF

!       DO jl = 1, nl
!         ! Longitude - OCLORG and FL are gridbbox centres
!         xjlon(jl) = (fl(jl)-oclorg+ocscal*0.5_dp)/ocscal+1
!         jlon = INT(xjlon(jl)+0.00001_dp)
!         IF (jlon <= 0) THEN
!            xjlon(jl) = xjlon(jl)+nlon
!            jlon = INT(xjlon(jl)+0.00001_dp)
!         ENDIF
!         jjlon(jl) = jlon
!       END DO

!       ! laction: Flag of points, where fdat still has to be computed
! !CDIR NODEP
!       DO jl = 1, nl
!         jlat = jjlat(jl)
!         jlon = jjlon(jl)
!         ! HD Land but OA Water
!         laction(jl) = (flag(jl,jb) > 0.5_dp .AND. foclsm(jlon,jlat) < 0.5_dp)
!       END DO

!       DO jl = 1, nl
!         jlat = jjlat(jl)
!         jlon = jjlon(jl)
!         ! HD Land but OA Water
!         IF (laction(jl)) THEN
          
!           ! Considered neighbour gridboxes in OA grid
!           ! N,S,W,E-Directions
!           xx1(jl) = 0.0_dp
!           xx2(jl) = 0.0_dp
!           IF (jlon /= 1) THEN
!             IF (foclsm(jlon-1,jlat) > 0.5_dp) THEN
!               xx1(jl) = fatmos(jlon-1,jlat)
!               xx2(jl) = 1.0_dp
!             ENDIF
!           ELSE
!             IF (foclsm(nlon,jlat) > 0.5_dp) THEN
!               xx1(jl) = fatmos(nlon,jlat)
!               xx2(jl) = 1.0_dp
!             ENDIF
!           ENDIF
!           IF (jlon /= nlon) THEN
!             IF (foclsm(jlon+1,jlat) > 0.5_dp) THEN
!               xx1(jl) = xx1(jl)+fatmos(jlon+1,jlat)
!               xx2(jl) = xx2(jl)+1.0_dp
!             ENDIF
!           ELSE
!             IF (foclsm(1,jlat) > 0.5_dp) THEN
!               xx1(jl) = xx1(jl)+fatmos(1,jlat)
!               xx2(jl) = xx2(jl)+1.0_dp
!             ENDIF
!           ENDIF
!           IF (jlat /= 1) THEN
!             IF (foclsm(jlon,jlat-1) > 0.5_dp) THEN
!               xx1(jl) = xx1(jl)+fatmos(jlon,jlat-1)
!               xx2(jl) = xx2(jl)+1.0_dp
!             ENDIF
!           ENDIF
!           IF (jlat /= nlat) THEN
!             IF (foclsm(jlon,jlat+1) > 0.5_dp) THEN
!               xx1(jl) = xx1(jl)+fatmos(jlon,jlat+1)
!               xx2(jl) = xx2(jl)+1.0_dp
!             ENDIF
!           ENDIF
!           ! Land point found?
!           IF (xx2(jl) > 0.5_dp) THEN 
!             fdat(jl,jb) = xx1(jl)/xx2(jl)
!             laction(jl) = .FALSE.
!           END IF
!         END IF
!       END DO
      
!       DO jl = 1, nl
!         jlat = jjlat(jl)
!         jlon = jjlon(jl)
!         ! HD Land but OA Water
!         IF (laction(jl)) THEN
!           ! if not --> NW,NE,SW,SE-Directions
!           IF (jlon /= 1) THEN
!             IF (jlat /= 1) THEN
!               IF (foclsm(jlon-1,jlat-1) > 0.5_dp) THEN
!                 xx1(jl) = xx1(jl)+fatmos(jlon-1,jlat-1)
!                 xx2(jl) = xx2(jl)+1.0_dp
!               ENDIF
!             ENDIF
!             IF (jlat /= nlat) THEN
!               IF (foclsm(jlon-1,jlat+1) > 0.5_dp) THEN
!                 xx1(jl) = xx1(jl)+fatmos(jlon-1,jlat+1)
!                 xx2(jl) = xx2(jl)+1.0_dp
!               ENDIF
!             ENDIF
!           ELSE
!             IF (jlat /= 1) THEN
!               IF (foclsm(nlon,jlat-1) > 0.5_dp) THEN
!                 xx1(jl) = xx1(jl)+fatmos(nlon,jlat-1)
!                 xx2(jl) = xx2(jl)+1.0_dp
!               ENDIF
!             ENDIF
!             IF (jlat /= nlat) THEN
!               IF (foclsm(nlon,jlat+1) > 0.5_dp) THEN
!                 xx1(jl) = xx1(jl)+fatmos(nlon,jlat+1)
!                 xx2(jl) = xx2(jl)+1.0_dp
!               ENDIF
!             ENDIF
!           ENDIF
!           IF (jlon /= nlon) THEN
!             IF (jlat /= 1) THEN
!               IF (foclsm(jlon+1,jlat-1) > 0.5_dp) THEN
!                 xx1(jl) = xx1(jl)+fatmos(jlon+1,jlat-1)
!                 xx2(jl) = xx2(jl)+1.0_dp
!               ENDIF
!             ENDIF
!             IF (jlat /= nlat) THEN
!               IF (foclsm(jlon+1,jlat+1) > 0.5_dp) THEN
!                 xx1(jl) = xx1(jl)+fatmos(jlon+1,jlat+1)
!                 xx2(jl) = xx2(jl)+1.0_dp
!               ENDIF
!             ENDIF
!           ELSE
!             IF (jlat /= 1) THEN
!               IF (foclsm(1,jlat-1) > 0.5_dp) THEN
!                 xx1(jl) = xx1(jl)+fatmos(1,jlat-1)
!                 xx2(jl) = xx2(jl)+1.0_dp
!               ENDIF
!             ENDIF
!             IF (jlat /= nlat) THEN
!               IF (foclsm(1,jlat+1) > 0.5_dp) THEN
!                 xx1(jl) = xx1(jl)+fatmos(1,jlat+1)
!                 xx2(jl) = xx2(jl)+1.0_dp
!               ENDIF
!             ENDIF
!           ENDIF
!           IF (xx2(jl) > 0.5_dp) THEN
!             ! Second next points in OA grid in N,S,W,E-Directions
!             fdat(jl,jb) = xx1(jl)/xx2(jl)
!             laction(jl) = .FALSE.
!           END IF
!         END IF
!       END DO
      
!       IF( ANY(laction) )   THEN
        
!         DO idd = 2, ndd
!           DO jl = 1, nl
!             jlat = jjlat(jl)
!             jlon = jjlon(jl)
!             ! HD Land but OA Water
!             IF (laction(jl)) THEN
              
!               xx1(jl) = 0.0_dp
!               xx2(jl) = 0.0_dp
!               IF (jlon-idd >= 1) THEN
!                 IF (foclsm(jlon-idd,jlat) > 0.5_dp) THEN
!                   xx1(jl) = fatmos(jlon-idd,jlat)
!                   xx2(jl) = 1.0_dp
!                 ENDIF
!               ELSE
!                 IF (foclsm(nlon+jlon-idd,jlat) > 0.5_dp) THEN
!                   xx1(jl) = fatmos(nlon+jlon-idd,jlat)
!                   xx2(jl) = 1.0_dp
!                 ENDIF
!               ENDIF
!               IF (jlon+idd <= nlon) THEN
!                 IF (foclsm(jlon+idd,jlat) > 0.5_dp) THEN
!                   xx1(jl) = xx1(jl)+fatmos(jlon+idd,jlat)
!                   xx2(jl) = xx2(jl)+1.0_dp
!                 ENDIF
!               ELSE
!                 IF (foclsm(jlon+idd-nlon,jlat) > 0.5_dp) THEN
!                   xx1(jl) = xx1(jl)+fatmos(jlon+idd-nlon,jlat)
!                   xx2(jl) = xx2(jl)+1.0_dp
!                 ENDIF
!               ENDIF
!               IF (jlat-idd >= 1) THEN
!                 IF (foclsm(jlon,jlat-idd) > 0.5_dp) THEN
!                   xx1(jl) = xx1(jl)+fatmos(jlon,jlat-idd)
!                   xx2(jl) = xx2(jl)+1.0_dp
!                 ENDIF
!               ENDIF
!               IF (jlat+idd <= nlat) THEN
!                 IF (foclsm(jlon,jlat+idd) > 0.5_dp) THEN
!                   xx1(jl) = xx1(jl)+fatmos(jlon,jlat+idd)
!                   xx2(jl) = xx2(jl)+1.0_dp
!                 ENDIF
!               ENDIF
!               ! End of Do (IDD) -Loop for Land Point found
!               IF (xx2(jl) > 0.5_dp) THEN
!                 fdat(jl,jb) = xx1(jl)/xx2(jl)
!                 laction(jl) = .FALSE.
!               END IF
!             END IF
!           ENDDO
          
!           ! end of HD land, OA water
!         ENDDO
!       END IF
!       IF( ANY(laction) )   THEN
!         DO jl = 1, nl
!           IF (laction(jl)) THEN
!             IF (ique /= 0) THEN
!               WRITE(message_text,*) 'no land point found for jl=',jl,  &
!                    '  jb=',jb, ': fl=',fl(jl), '  fb=',fb(jl)
!               CALL message('hydrology_corr', message_text)
!             END IF
!           END IF
!         END DO
!       END IF
!     ENDDO
    
!     x1 = 0.0_dp
!     x2 = 0.0_dp
!     DO jl = 1, nlon
!       x1 = x1+SUM(fatmos(jl,:)*foclsm(jl,:)*aoarea(:))
!     ENDDO
!     DO jl = 1, nl
!       x2 = x2+SUM(fdat(jl,:)*flag(jl,:)*area(:))
!     ENDDO
!     xresi = xresi+x1-x2
    

!   END SUBROUTINE hydrology_corr
  SUBROUTINE hydrology_corr(nlon, nlat, fatmos, foclsm, aoarea, philat, &
                            fdat, hd_lsm, area, xresi,                   &
                            oclorg, ocscal, florg, fborg, fscal, lhd_que)


    ! **** Correction of the atmospheric grid -> 0.5 degree transformation
    !
    ! ***** Version 1.0 - November 1999
    ! Programmed and developed by Stefan Hagemann, MPI
    ! Remark: Input data FATMOS (Runoff or Drainage) should have the unit m/s.
    !
    ! ***** Version 1.1 - January 2001
    ! Longitude-Index-Correction
    !
    ! ***** Version 2.1 - January 2001
    ! ECHAM5- Version incl. Gaussian latitudes
    !
    !
    ! **** Global Arrays:
    !
    !  fdat = Data array at 0.5 degree
    !  hd_lsm = Land mask array at 0.5 degree
    ! area(jb) = Array of Gridbox areas, Unit = [m^2]
    !
    ! ***** Parameter and arrays for atmosphere ocean grid
    !
    !    nlon = Longitudes of global ocean grid
    !    nlat = Latitudes of global ocean grid
    !
    !  oclorg = Longitudinal origin of global ocean grid
    !  ocscal = resolution/scale = Width of an ocean Gridbox in degree
    !
    !  fatmos = atmospheric data array
    !  foclsm = Land Sea Mask on Ocean grid
    !  aoarea = Area per gridbox per latitude of atmospheric grid [m^2]
    !   philat = Gaussian latitude of global ocean grid (centre coordinates)
    !
    !   xresi = Residuum (Runoff+Drainage), which results from different
    !           land sea masks of the atmospheric grid and the 0.5 degree grid.
    !           In the latest HD model version, a start value is passed to the
    !           HD model that may include further residual water terms that should
    !           distributed with the discharge to close the water balance in the
    !           coupled atmosphere ocean system.
    !
    !    lhd_que = Log-Output switch ( 0 = No Log-output to STDOUT)
    !
    INTEGER, INTENT(in) :: nlon, nlat
    ! Input fields AO-Grid

    REAL(dp), INTENT(in) :: foclsm(nlon,nlat), aoarea(nlat)
    REAL(dp), INTENT(in) :: fatmos(nlon,nlat), philat(nlat)
    REAL(dp), INTENT(inout) :: fdat(nl,nb)
    REAL(sp), INTENT(in) :: hd_lsm(nl,nb)
    REAL(dp), INTENT(in) :: area(nb)
    REAL(dp), INTENT(inout) :: xresi
    REAL(dp), INTENT(in) :: oclorg, ocscal
    REAL(dp), INTENT(in) :: florg, fborg, fscal
    LOGICAL, INTENT(in) :: lhd_que

    REAL(dp) :: x1, x2
    INTEGER :: jb, jl, jlon, jlat

    REAL(dp) :: fb, fl(nl)
    INTEGER :: jjlon(nl)
    LOGICAL :: reassigned


    ! OA = Ocean-Atmosphere Grid, HD = 0.5 Grad HD-Model Grid
    DO jl = 1, nl
      fl(jl) = MOD(REAL(jl, dp) * fscal + florg - 0.5_dp * fscal + 180._dp, &
           fullcirc) - 180._dp
    END DO

    DO jl = 1, nl
      ! Longitude - OCLORG and FL are gridbbox centres
      jjlon(jl) = INT((MOD(fl(jl) - oclorg + ocscal*0.5_dp - fullcirc, &
           -fullcirc) + fullcirc)/ocscal + 1 + 0.00001_dp)
    END DO

    DO jb = 1, nb-6
      ! First, compute fb and fl for all jl values
      ! grid box centre:
      fb = fborg-REAL(jb,dp)*fscal+0.5_dp*fscal
      ! Corresponding Index in Ocean Grid
      ! Latitude

!!! without Gauss:    XJLAT = (OCBORG-FB) / OCSCAL+1.  , OCBORG = Borderline
!!!                 JLAT = INT(XJLAT+0.00001)

      ! and compute jlat and jlog for all jl values

      ! mz_bk_20120730+
      ! jlat = dec_monotonic_interval_closest_midpoint(philat, &
      !      fb, aub=90._dp, alb=-90._dp)
      ! mz_bk_20120730-
      jlat = dec_mon_interval_clos_midpoint(philat, &
           fb, aub=90._dp, alb=-90._dp)
      ! laction: Hd_lsm of points, where fdat still has to be computed
      DO jl = 1, nl
        ! HD Land but OA Water?
        IF (hd_lsm(jl,jb) > 0.5_sp &
             .AND. foclsm(jjlon(jl), jlat) < 0.5_dp) THEN
          jlon = jjlon(jl)
          ! if not --> NW,NE,SW,SE-Directions
          reassigned = reassign_runoff(nlon, nlat, jlon, jlat, &
               foclsm, fatmos, fdat(jl, jb))
          IF (.NOT. reassigned .AND. lhd_que) THEN
            WRITE(message_text,*) 'no land point found for jl=',jl,  &
                 '  jb=',jb, ': fl=',fl(jl), '  fb=',fb
            ! mz_bk_20120730+
            ! CALL message('hydrology_corr', message_text)
            ! mz_bk_20120730-
          END IF
        END IF
      END DO
    ENDDO

    x1 = 0.0_dp
    x2 = 0.0_dp
    DO jl = 1, nlon
      x1 = x1+SUM(fatmos(jl,:)*foclsm(jl,:)*aoarea(:))
    ENDDO
    DO jl = 1, nl
      x2 = x2+SUM(fdat(jl,:) * REAL(hd_lsm(jl,:), dp) * area(:))
    ENDDO
    xresi = xresi+x1-x2

  END SUBROUTINE hydrology_corr

  FUNCTION reassign_runoff(nlon, nlat, jlon, jlat, foclsm, fatmos, fdat) &
       RESULT(reassigned)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nlon, nlat, jlon, jlat
    REAL(dp), INTENT(in) :: foclsm(nlon, nlat), fatmos(nlon, nlat)
    REAL(dp), INTENT(inout) :: fdat
    LOGICAL :: reassigned

    REAL(dp) :: x1, inverted_neighbour_weight
    INTEGER :: idd

    reassigned = .FALSE.
    ! HD Land but OA Water
    ! Considered neighbour gridboxes in OA grid
    ! N,S,W,E-Directions
    x1 = 0.0_dp
    inverted_neighbour_weight = 0.0_dp
    IF (jlon /= 1) THEN
      IF (foclsm(jlon-1,jlat) > 0.5_dp) THEN
        x1 = fatmos(jlon-1,jlat)
        inverted_neighbour_weight = 1.0_dp
      ENDIF
    ELSE
      IF (foclsm(nlon,jlat) > 0.5_dp) THEN
        x1 = fatmos(nlon,jlat)
        inverted_neighbour_weight = 1.0_dp
      ENDIF
    ENDIF
    IF (jlon /= nlon) THEN
      IF (foclsm(jlon+1,jlat) > 0.5_dp) THEN
        x1 = x1+fatmos(jlon+1,jlat)
        inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
      ENDIF
    ELSE
      IF (foclsm(1,jlat) > 0.5_dp) THEN
        x1 = x1+fatmos(1,jlat)
        inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
      ENDIF
    ENDIF
    IF (jlat /= 1) THEN
      IF (foclsm(jlon,jlat-1) > 0.5_dp) THEN
        x1 = x1+fatmos(jlon,jlat-1)
        inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
      ENDIF
    ENDIF
    IF (jlat /= nlat) THEN
      IF (foclsm(jlon,jlat+1) > 0.5_dp) THEN
        x1 = x1+fatmos(jlon,jlat+1)
        inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
      ENDIF
    ENDIF
    ! Land point found?
    IF (inverted_neighbour_weight > 0.5_dp) THEN
      fdat = x1/inverted_neighbour_weight
      reassigned = .TRUE.
      RETURN
    END IF

    IF (jlon /= 1) THEN
      IF (jlat /= 1) THEN
        IF (foclsm(jlon-1,jlat-1) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon-1,jlat-1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(jlon-1,jlat+1) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon-1,jlat+1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
    ELSE
      IF (jlat /= 1) THEN
        IF (foclsm(nlon,jlat-1) > 0.5_dp) THEN
          x1 = x1+fatmos(nlon,jlat-1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(nlon,jlat+1) > 0.5_dp) THEN
          x1 = x1+fatmos(nlon,jlat+1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
    ENDIF
    IF (jlon /= nlon) THEN
      IF (jlat /= 1) THEN
        IF (foclsm(jlon+1,jlat-1) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon+1,jlat-1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(jlon+1,jlat+1) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon+1,jlat+1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
    ELSE
      IF (jlat /= 1) THEN
        IF (foclsm(1,jlat-1) > 0.5_dp) THEN
          x1 = x1+fatmos(1,jlat-1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(1,jlat+1) > 0.5_dp) THEN
          x1 = x1+fatmos(1,jlat+1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
    ENDIF
    IF (inverted_neighbour_weight > 0.5_dp) THEN
      ! Second next points in OA grid in N,S,W,E-Directions
      fdat = x1/inverted_neighbour_weight
      reassigned = .TRUE.
      RETURN
    END IF
    extended_surround_loop: DO idd = 2, ndd
      ! HD Land but OA Water
      x1 = 0.0_dp
      inverted_neighbour_weight = 0.0_dp
      IF (jlon-idd >= 1) THEN
        IF (foclsm(jlon-idd,jlat) > 0.5_dp) THEN
          x1 = fatmos(jlon-idd,jlat)
          inverted_neighbour_weight = 1.0_dp
        ENDIF
      ELSE
        IF (foclsm(nlon+jlon-idd,jlat) > 0.5_dp) THEN
          x1 = fatmos(nlon+jlon-idd,jlat)
          inverted_neighbour_weight = 1.0_dp
        ENDIF
      ENDIF
      IF (jlon+idd <= nlon) THEN
        IF (foclsm(jlon+idd,jlat) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon+idd,jlat)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ELSE
        IF (foclsm(jlon+idd-nlon,jlat) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon+idd-nlon,jlat)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat-idd >= 1) THEN
        IF (foclsm(jlon,jlat-idd) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon,jlat-idd)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat+idd <= nlat) THEN
        IF (foclsm(jlon,jlat+idd) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon,jlat+idd)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      ! End of Do (IDD) -Loop for Land Point found
      IF (inverted_neighbour_weight > 0.5_dp) THEN
        fdat = x1/inverted_neighbour_weight
        reassigned = .TRUE.
        EXIT extended_surround_loop
      END IF
    END DO extended_surround_loop
    ! end of HD land, OA water
  END FUNCTION reassign_runoff

  ! mz_bk_20120730-


!-------------------------------------------------------------------------------

  ! mz_bk_20120730+
  ! SUBROUTINE kasglob(finp, ymod, a_k, a_n, iflow, fmem, mm)
  !   !
  !   ! ***** Global Flow Simulation with the conceptual model reservoir cascade
  !   !   Program was partailly written following the routine lfsim (in gate.for)
  !   !   and funkas (in modfunct.for).
  !   !
  !   ! ***** Programmed and developed by Stefan Hagemann
  !   !
  !   ! ***** Version 2.2 - November 1995
  !   !   Finer resolution of system function for Riverflow
  !   !   Re-arranging of IN/OUTput-arrays and Passing within a 
  !   !      single reservoir field fmem
  !   !
  !   ! ***** Version 3.0 - November 1995
  !   !   Computation of outflow via Differential Equation
  !   !   of lin. reservoir cascade, comprising of nn reservoirs with
  !   !   nn = INT(a_n) = INT(n)
  !   !
  !   ! ***** Version 3.1 - Februar 1996
  !   !   Implementation of possible computation of riverflow with mm
  !   !   sub (internal) time steps.
  !   !
  !   ! ***** Version 5.0 - Oktober 1999
  !   !   Runoff-Input-data are passed to kasglob instead of reading it within 
  !   !   kasglob itself.
  !   !   Calling parameters/variables luinp and area are deleted.
  !   !
  !   !
  !   ! ***** List of Variables:
  !   !
  !   !    ii = Number of time step
  !   !  finp = Input array for Overlandflow / Riverflow
  !   !  ymod = simulated Overlandflow / Riverflow Array
  !   !   a_k = Array of  k-Parameter [day]
  !   !   a_n = Array of n-Parameter
  !   ! iflow = Flow type
  !   !         1 = Overlandflow
  !   !         2 = Riverflow
  !   !
  !   !  fmem(nl, nb, nmem) = Intermediate content of reservoir cascade
  !   !                       for Inflows or Runoffs per Gridbox
  !   !
  !   !    mm = Computation of riverflow in mm Sub time steps
  !   !
  !   ! ***** Local Variables
  !   !
  !   !  fdum = Dummy
  !   ! akdiv = Dummy factor for division
  !   !
  !   !

  !   INTEGER, INTENT(in) :: iflow, mm
  !   REAL(dp), INTENT(in) :: a_k(nl, nb), a_n(nl, nb),  finp(nl, nb)
  !   REAL(dp), INTENT(inout) :: ymod(nl, nb)
  !   REAL(dp), INTENT(inout) :: fmem(:,:,:)

  !   REAL(dp) :: akdiv, fdum, amod, divmm, fakmm

  !   INTEGER :: j, jl, jb
  !   INTEGER :: nn

  !   ! Transformation of Input values of time step II and
  !   ! temporary storage in array FINP

  !   IF (iflow == 1) THEN
  !     divmm = 1.0_dp
  !     fakmm = 1.0_dp
  !   ELSE IF (iflow == 2) THEN
  !     divmm = 1.0_dp/REAL(mm,dp)
  !     fakmm = REAL(mm,dp)
  !   ENDIF


  !   ! **** Computing modeled value at each grid point

  !   DO jb = 1, nb
  !     DO jl = 1, nl
  !       IF (a_n(jl,jb) > 0.5_dp) THEN

  !         nn = INT(a_n(jl,jb)+0.5_dp)
  !         amod = a_k(jl,jb)*a_n(jl,jb)/REAL(nn,dp)

  !         ! *** Dt=1 day ==> AKDIV = 1./ (AMOD+1.)
  !         akdiv = 1.0_dp/(amod+divmm)

  !         ! *** Input at time step II
  !         fdum = finp(jl,jb)*divmm

  !         ! *** Nash-Cascade
  !         ! *** It is [AMOD] = day ==> AMOD(sec) = AMOD(day) * 1 day
  !         ! *** Remember: In principle,it is: FDUM = FINP * 1 day
  !         ! ***           ==> FMEM = x * 1 day
  !         ! ***           ==> FDUM = x * 1 day * 1 / AMOD(sec)
  !         ! ***                    = x * 1 day / (AMOD(day) * 1 day)
  !         ! ***                    = x / AMOD(day)
  !         ! ***           ==> FMEM = x * 1 day - FDUM * 1 day
  !         ! ***                    = (x - FDUM) * 1 day
  !         ! *** Outflow FDUM is computed correctly, Intermediate reservoir unit is
  !         ! *** a volume flow instead of a volume. This is to avoid 
  !         ! *** back and forth multiplication with factor 1 day = 86400 sec

  !         DO j = 1, nn
  !           fmem(jl,jb,j) = fmem(jl,jb,j)+fdum
  !           fdum = fmem(jl,jb,j)*akdiv*divmm
  !           fmem(jl,jb,j) = fmem(jl,jb,j)-fdum
  !         ENDDO
  !         ymod(jl,jb) = fdum*fakmm

  !       ENDIF
  !     ENDDO
  !   ENDDO

  ! END SUBROUTINE kasglob

  SUBROUTINE kasglob(finp, ymod, a_k, a_n, iflow, fmem, a_n_kas)
    !
    ! ***** Global Flow Simulation with the conceptual model reservoir cascade
    !   Program was partailly written following the routine lfsim (in gate.for)
    !   and funkas (in modfunct.for).
    !
    ! ***** Programmed and developed by Stefan Hagemann
    !
    ! ***** Version 2.2 - November 1995
    !   Finer resolution of system function for Riverflow
    !   Re-arranging of IN/OUTput-arrays and Passing within a
    !      single reservoir field fmem
    !
    ! ***** Version 3.0 - November 1995
    !   Computation of outflow via Differential Equation
    !   of lin. reservoir cascade, comprising of nn reservoirs with
    !   nn = INT(a_n) = INT(n)
    !
    ! ***** Version 3.1 - Februar 1996
    !   Implementation of possible computation of riverflow with mm
    !   sub (internal) time steps.
    !
    ! ***** Version 5.0 - Oktober 1999
    !   Runoff-Input-data are passed to kasglob instead of reading it within
    !   kasglob itself.
    !   Calling parameters/variables luinp and area are deleted.
    !
    !
    ! ***** List of Variables:
    !
    !    ii = Number of time step
    !  finp = Input array for Overlandflow / Riverflow
    !  ymod = simulated Overlandflow / Riverflow Array
    !   a_k = Array of  k-Parameter [day]
    !   a_n = Array of n-Parameter
    ! iflow = Flow type
    !         1 = Overlandflow
    !         2 = Riverflow
    !
    !  fmem(nl, nb, nmem) = Intermediate content of reservoir cascade
    !                       for Inflows or Runoffs per Gridbox
    !
    ! ***** Local Variables
    !
    !  fdum = Dummy
    ! akdiv = Dummy factor for division
    !
    !

    INTEGER, INTENT(in) :: iflow
    REAL(dp), INTENT(in) :: a_k(nl, nb), a_n(nl, nb),  finp(nl, nb)
    REAL(dp), INTENT(inout) :: ymod(nl, nb)
    REAL(dp), INTENT(inout) :: fmem(:,:,:)
    TYPE(cart_xidx_2d), INTENT(in) :: a_n_kas(:)

    REAL(dp) :: akdiv, fdum, amod, divmm, fakmm, fmd_sum

    INTEGER :: j, jl, jb
    INTEGER :: nn
    INTEGER :: nx, i, extlen, extelem

    ! Transformation of Input values of time step II and
    ! temporary storage in array FINP

    IF (iflow == overlandflow) THEN
      divmm = 1.0_dp
      fakmm = 1.0_dp
    ELSE IF (iflow == riverflow) THEN
      divmm = 1.0_dp/REAL(mm,dp)
      fakmm = REAL(mm,dp)
    ENDIF

    nx = SIZE(a_n_kas)
    ! **** Computing modeled value at each grid point

    DO i = 1, nx
      jb = a_n_kas(i)%ilat
      jl = a_n_kas(i)%ilon
      extlen = a_n_kas(i)%extlen
      DO extelem = 1, extlen
        nn = NINT(a_n(jl,jb))
        amod = a_n_kas(i)%amod(extelem)

        ! *** Dt=1 day ==> AKDIV = 1./ (AMOD+1.)
        akdiv = a_n_kas(i)%akdiv(extelem)

        ! *** Input at time step II
        fdum = finp(jl,jb)*divmm

        ! *** Nash-Cascade
        ! *** It is [AMOD] = day ==> AMOD(sec) = AMOD(day) * 1 day
        ! *** Remember: In principle,it is: FDUM = FINP * 1 day
        ! ***           ==> FMEM = x * 1 day
        ! ***           ==> FDUM = x * 1 day * 1 / AMOD(sec)
        ! ***                    = x * 1 day / (AMOD(day) * 1 day)
        ! ***                    = x / AMOD(day)
        ! ***           ==> FMEM = x * 1 day - FDUM * 1 day
        ! ***                    = (x - FDUM) * 1 day
        ! *** Outflow FDUM is computed correctly, Intermediate reservoir unit is
        ! *** a volume flow instead of a volume. This is to avoid
        ! *** back and forth multiplication with factor 1 day = 86400 sec

        DO j = 1, nn
          fmd_sum = fmem(jl,jb,j) + fdum
          fdum = fmd_sum * akdiv * divmm
          fmem(jl,jb,j) = fmd_sum - fdum
        END DO
        ymod(jl,jb) = fdum*fakmm
        jl = jl + 1
      ENDDO
    ENDDO

  END SUBROUTINE kasglob

  ! mz_bk_20120730-

!-------------------------------------------------------------------------------

  ! mz_bk_20120730+
!   SUBROUTINE hydrology_to_ocean(nlon, nlat, oclorg, ocscal, fglat,  &
!                                 friv, fdir, disch, foclsm, xresi, ique )

!     !
!     ! ******* This programs distributes the river discharge from the 0.5
!     !         degree inflow points into the ocean into the considered
!     !         ocean gridbox
!     !
!     ! **** originally fixed for Ocean-Grid=T42
!     !
!     !  Programmed and developed by Stefan Hagemann, MPI
!     !
!     ! ***** Version 1.0 -- Oktober 1999
!     !
!     ! ***** Version 2.0 -- January 2001
!     !     ECHAM5- Version incl. Gaussian latitudes
!     !
!     ! ****** List of Variables
!     !
!     !  friv = Inflow array on HD model grid
!     !  fdir = River direction file that defines river mouthes (destinations) as 0
!     !         on HD model grid
!     !  xidb = Summation array of inflows, for which no inflowbox into the 
!     !         ocean was found, e.g. Kaspian Sea and Interior
!     !         Drainage Basins
!     !  kfound = Inflow-Point found on ocean grid YES/NO
!     !
!     !  nlon = Longitudes of global ocean grid
!     !  nlat = Latitudes of global ocean grid
!     !  oclorg = Longitudinal origin of global ocean grid
!     !  ocscal = Scale/Resolution = Width of Ocean Gridbox in degree
!     !   fglat = Gaussian latitude of global ocean grid (centre coordinates)
!     !
!     !  foclsm = Land Sea Mask on Ocean grid
!     !  kfocim = Mask of inflow points on Ocean grid
!     !   disch = Inflow array on Ocean grid
!     !   xresi = Residuum (Runoff+Drainage), which results from different 
!     !           land sea masks of the atmospheric grid and the 0.5 degree grid.
!     !           In the latest HD model version, a start value is passed to the
!     !           HD model that may include further residual water terms that should
!     !           distributed with the discharge to close the water balance in the
!     !           coupled atmosphere ocean system.
!     !           XIDB is added to Xresi.
!     !    ique = Log-Output switch ( 0 = No Log-output to STDOUT)
!     !

!     INTEGER, INTENT(in) :: nlon, nlat

!     REAL(dp) :: friv(nl,nb), fdir(nl,nb)
!     REAL(dp) :: disch(nlon,nlat), foclsm(nlon,nlat), fglat(nlat)
!     REAL(dp) :: oclorg, ocscal, xresi
!     INTEGER :: kfocim(nlon,nlat)
!     REAL(dp) :: xanf, xend, xjlon, xjlat, xidb, fb, fl
!     INTEGER :: ique,jl,jb, kfound, jlat, jlon

!     disch(:,:) = 0.0_dp
!     kfocim(:,:) = 0
!     xidb = xresi

!     ! ******* Loop over all inflow points

!     DO jb = 1, nb
!       DO jl = 1, nl
!         IF (fdir(jl,jb) == 0 .OR. fdir(jl,jb) == 5) THEN

!           ! *** Gridbox centre 

!           fb = fborg-REAL(jb,dp)*fscal+0.5_dp*fscal
!           fl = REAL(jl,dp)*fscal+florg-0.5_dp*fscal
!           IF (fl >= 180.0_dp) fl = fl-360.0_dp

!           ! *** Corrwsponding Index in Ocean Grid
!           !
!           ! *** Latitude
! !!!  without Gauss: XJLAT = (OCBORG-FB)/OCSCAL+1, OCBORG = Borderline
! !!!               JLAT = INT(XJLAT+0.00001)
!           !
!           DO jlat = 1, nlat
!             IF (jlat == 1) THEN
!               xanf = 90.0_dp
!             ELSE
!               xanf = (fglat(jlat-1)+fglat(jlat))*0.5_dp
!             ENDIF
!             IF (jlat == nlat) THEN
!               xend = -90.0_dp
!             ELSE
!               xend = (fglat(jlat+1)+fglat(jlat))*0.5_dp
!             ENDIF
!             IF (fb <= xanf .AND. fb > xend) THEN
!               xjlat = jlat+(xanf-fb)/(xanf-xend)
!               EXIT
!             ENDIF
!           ENDDO
!           IF (jlat == nlat+1) THEN
!             WRITE(message_text,*) ' error in jlat=', jlat
!             CALL message('hydrology_to_ocean', message_text)
!             jlat = nlat
!             xjlat = jlat+(xanf-fb)/(xanf-xend)
!           ENDIF

!           ! *** Longitude - OCLORG and FL are Gridbox centres

!           xjlon = (fl-oclorg+ocscal*0.5_dp)/ocscal+1.0_dp
!           jlon = INT(xjlon+0.00001_dp)
!           IF (jlon <= 0) THEN
!             xjlon = xjlon+nlon
!             jlon = INT(xjlon+0.00001_dp)
!           ENDIF

!           ! *** Mouth/destination point = Ocean Point in Ocean grid?
!           IF (fdir(jl,jb) == 0) THEN
!             IF (foclsm(jlon,jlat) < 0.5_dp) THEN
!               disch(jlon,jlat) = disch(jlon,jlat)+friv(jl,jb)
!               kfocim(jlon,jlat) = 1

!               ! *** searching for closest Ocean Point in Ocean grid
!             ELSE
!               IF (ique /= 0)  THEN
!                 WRITE(message_text,*) 'jlon = ', jlon, '  jlat = ', jlat
!                 CALL message('hydrology_to_ocean', message_text)                
!               END IF
!               CALL oclook(foclsm, nlon,nlat, jlon,jlat, xjlon, xjlat, kfound)
!               IF (kfound == 1) THEN
!                 disch(jlon,jlat) = disch(jlon,jlat)+friv(jl,jb)
!                 kfocim(jlon,jlat) = 1
!                 IF (ique /= 0) THEN
!                   WRITE(message_text,*) '--> jlon = ', jlon, '  jlat = ', jlat
!                   CALL message('hydrology_to_ocean', message_text)
!                 END IF
!               ELSE
!                 IF (ique /= 0) THEN
!                   WRITE(message_text,*) 'no ocean point at (jl,jb)=',         &
!                        jl,jb, ' --> jlon = ', jlon, '  jlat = ', jlat
!                   CALL message('hydrology_to_ocean', message_text)
!                 END IF
!                 xidb = xidb+friv(jl,jb)
!               ENDIF
!             ENDIF

!             !    *** interior drainage point?
!           ELSE IF (fdir(jl,jb) == 5) THEN
!             IF (ique /= 0) THEN
!               WRITE(message_text,*) 'idb at (jl,jb)=',                    &
!                    jl,jb, ' --> jlon = ', jlon, '  jlat = ', jlat
!               CALL message('hydrology_to_ocean', message_text)
!             END IF
!             xidb = xidb+friv(jl,jb)
!           ENDIF
!         ENDIF
!       ENDDO
!     ENDDO

!     ! Distributing the water in XIDB to all Ocean Inflow Points
!     ! Applying a weight to treat arid and humid regions differently

!     kfound = SUM(kfocim(:,:))

!     IF (kfound > 0) THEN
!       disch(:,:) = disch(:,:)+disch(:,:)/SUM(disch(:,:))*  &
!            xidb*REAL(kfocim(:,:),dp)
!     ELSE
!       WRITE(message_text,*) 'error no inflow points on ocean grid found'
!       CALL message('hydrology_to_ocean', message_text)

!     ENDIF
!   END SUBROUTINE hydrology_to_ocean
  SUBROUTINE hydrology_to_ocean(nlon, nlat, oclorg, ocscal, philat,  &
                                friv, fdir, disch, xresi, lhd_que )

    !
    ! ******* This programs distributes the river discharge from the 0.5
    !         degree inflow points into the ocean into the considered
    !         ocean gridbox
    !
    ! **** originally fixed for Ocean-Grid=T42
    !
    !  Programmed and developed by Stefan Hagemann, MPI
    !
    ! ***** Version 1.0 -- Oktober 1999
    !
    ! ***** Version 2.0 -- January 2001
    !     ECHAM5- Version incl. Gaussian latitudes
    !
    ! ****** List of Variables
    !
    !  friv = Inflow array on HD model grid
    !  fdir = River direction file that defines river mouthes (destinations) as 0
    !         on HD model grid
    !  xidb = Summation array of inflows, for which no inflowbox into the
    !         ocean was found, e.g. Kaspian Sea and Interior
    !         Drainage Basins
    !  any_ocinflow = Inflow-Point found on ocean grid .TRUE./.FALSE.
    !
    !  nlon = Longitudes of global ocean grid
    !  nlat = Latitudes of global ocean grid
    !  oclorg = Longitudinal origin of global ocean grid
    !  ocscal = Scale/Resolution = Width of Ocean Gridbox in degree
    !  philat = Gaussian latitude of global ocean grid (centre coordinates)
    !
    !   disch = Inflow array on Ocean grid
    !   xresi = Residuum (Runoff+Drainage), which results from different
    !           land sea masks of the atmospheric grid and the 0.5 degree grid.
    !           In the latest HD model version, a start value is passed to the
    !           HD model that may include further residual water terms that should
    !           distributed with the discharge to close the water balance in the
    !           coupled atmosphere ocean system.
    !           XIDB is added to Xresi.
    ! lhd_que = Log-Output switch (.FALSE. = No Log-output to STDOUT)
    !

    INTEGER, INTENT(in) :: nlon, nlat

    REAL(dp), INTENT(in) :: friv(nl,nb)
    INTEGER, INTENT(in) :: fdir(nl,nb)
    REAL(dp), INTENT(out) :: disch(nlon,nlat)
    REAL(dp), INTENT(in) :: philat(nlat)
    REAL(dp), INTENT(in) :: oclorg, ocscal, xresi
    LOGICAL, INTENT(in) :: lhd_que
#if 0
    REAL(dp) :: xjlat
#endif
    REAL(dp) :: xjlon, xidb, fb, fl
    INTEGER :: jl,jb, jlat, jlon
    TYPE(cart_idx_2d) :: dest
    LOGICAL :: any_ocinflow

    disch(:,:) = 0.0_dp
    xidb = xresi
    any_ocinflow = .FALSE.
    ! ******* Loop over all inflow points

    DO jb = 1, nb
      fb = fborg-REAL(jb,dp)*fscal+0.5_dp*fscal
      ! mz_bk_20120730+
      ! jlat = dec_monotonic_interval_closest_midpoint(philat, &
      !      fb, aub=90._dp, alb=-90._dp)
      ! mz_bk_20120730-
      jlat = dec_mon_interval_clos_midpoint(philat, &
           fb, aub=90._dp, alb=-90._dp)
      DO jl = 1, nl
        IF (fdir(jl,jb) == 0 .OR. fdir(jl,jb) == 5) THEN

          ! *** Gridbox centre

          fl = MOD(REAL(jl,dp)*fscal+florg-0.5_dp*fscal + 180._dp, &
               fullcirc) - 180.0_dp
          ! *** Corresponding Index in Ocean Grid
          !
          ! *** Latitude
!!!  without Gauss: XJLAT = (OCBORG-FB)/OCSCAL+1, OCBORG = Borderline
!!!               JLAT = INT(XJLAT+0.00001)
          !
!          xanf = MERGE((philat(jlat-1)+philat(jlat))*0.5_dp, 90.0_dp, jlat /= 1)
!          xend = MERGE((philat(jlat+1)+philat(jlat))*0.5_dp, -90.0_dp, &
!               jlat /= nlat)
!          xjlat = jlat + (xanf-fb)/(xanf-xend)
#if 0
          IF (jlat < 1 .OR. jlat > nlat) THEN
            WRITE(message_text,*) ' error in jlat=', jlat
            CALL message('hydrology_to_ocean', message_text)
            jlat = nlat
            xjlat = jlat+(xanf-fb)/(xanf-xend)
          ENDIF
#endif
          ! *** Longitude - OCLORG and FL are Gridbox centres

          xjlon = (fl - oclorg + ocscal * 0.5_dp)/ocscal + 1.0_dp
          jlon = INT(xjlon + 0.00001_dp)
          IF (jlon <= 0) THEN
            xjlon = xjlon + nlon
            jlon = INT(xjlon + 0.00001_dp)
          ENDIF

          ! *** Mouth/destination point = Ocean Point in Ocean grid?
          IF (fdir(jl, jb) == 0) THEN
            IF (gl_slf(jlon, jlat) + gl_alake(jlon, jlat) < 0.5_dp &
                 .AND. gl_slm(jlon, jlat) <= 0.5_dp) THEN
              disch(jlon, jlat) = disch(jlon, jlat) + friv(jl, jb)
              any_ocinflow = .TRUE.
              ! *** searching for closest Ocean Point in Ocean grid
            ELSE
              IF (lhd_que)  THEN
                WRITE(message_text,*) 'jlon = ', jlon, '  jlat = ', jlat
                CALL message('hydrology_to_ocean', message_text)
              END IF
              dest = oclook_cache(jl, jb)
              IF (dest%ilon /= -1) THEN
                disch(dest%ilon,dest%ilat) = disch(dest%ilon,dest%ilat) &
                  + friv(jl,jb)
                any_ocinflow = .TRUE.
                IF (lhd_que) THEN
                  WRITE(message_text,*) '--> jlon = ', dest%ilon, &
                       '  jlat = ', dest%ilat
                  CALL message('hydrology_to_ocean', message_text)
                END IF
              ELSE
                IF (lhd_que) THEN
                  WRITE(message_text,*) 'no ocean point at (jl,jb)=',         &
                       jl,jb, ' --> jlon = ', jlon, '  jlat = ', jlat
                  CALL message('hydrology_to_ocean', message_text)
                END IF
                xidb = xidb + friv(jl,jb)
              ENDIF
            ENDIF

            !    *** interior drainage point?
          ELSE ! IF (fdir(jl,jb) == 5) THEN
            IF (lhd_que) THEN
              WRITE(message_text,*) 'idb at (jl,jb)=',                    &
                   jl,jb, ' --> jlon = ', jlon, '  jlat = ', jlat
              CALL message('hydrology_to_ocean', message_text)
            END IF
            xidb = xidb + friv(jl,jb)
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    ! Distributing the water in XIDB to all Ocean Inflow Points
    ! Applying a weight to treat arid and humid regions differently

    IF (any_ocinflow) THEN
      disch(:,:) = disch(:,:)+disch(:,:)/SUM(disch(:,:)) * xidb
    ELSE
      WRITE(message_text,*) 'error no inflow points on ocean grid found'
      CALL message('hydrology_to_ocean', message_text)

    ENDIF
  END SUBROUTINE hydrology_to_ocean

  ! mz_bk_20120730-

!-------------------------------------------------------------------------------

  ! mz_bk_20120730+
!   SUBROUTINE glacier_to_ocean
    
!     USE messy_main_data_bi,      ONLY: ngl, nlon, lcouple,gridarea
!     USE messy_main_timer_bi,     ONLY: get_time_step => timer_get_time_step ! um_ak_20090625
!     USE messy_main_tools,        ONLY: find_next_free_unit

!     !   ***** Glacier calving subroutine fuer ECHAM5
!     !         Distributes the P-E differences over glacier grid boxes 
!     !         to the nearest ocean grid box indicated by the land sea mask slm.
!     !         For accurate distribution, the gridbox area of source and target
!     !         grid box have to be taken into account.
!     !
!     !    Programmed and developed by Stefan Hagemann, MPI
!     !
!     !   ***** Version 1.0 -- February 2001
!     !         Programmierung analogous to hd_tooc.f90
!     !         Subroutine uses oclook.f that is included in hd_tooc.f
!     !
!     !   ****** Variablenliste
!     !   
!     !     apmecal = P minus E array -- expected in [m/s]
!     !      xidb = Summation array of Inflows, for which no Inflowbox in 
!     !             Ocean grid was found.
!     !    kfound = Inflow-Point found  on ocean grid YES/NO
!     !   
!     !      nlon = Longitudes of global atmosphere/ocean grid
!     !       ngl = Latitudes of global ocean grid
!     !   
!     !       slm = Land Sea Mask on atmosphere/ocean grid
!     !    kfocim = Mask of inflow points on Ocean grid
    
!     LOGICAL :: lex
    
!     INTEGER :: kfocim(nlon,ngl)
!     INTEGER :: ique, jl, jb, jlat, jlon, kfound
!     INTEGER :: iunit, istep

! !    REAL(dp) :: zcalv(nlon, ngl)   
!     REAL(dp) :: xjlon, xjlat, xidb
! !    REAL(dp) :: zd, zarea

!     ique = 0

!     gl_zcalv(:,:)  = 0.0_dp
!     kfocim(:,:) = 0
!     xidb        = 0.0_dp

! !    zd    = 0.0_dp
! !    zarea = 0.0_dp

!     !    Redrisbution of P-E values over glacier points to closest
!     !    ocean points
!     !
!     !   ******* Loop over all glacier points

!     DO jb= 1, ngl
!       DO jl= 1, nlon
!         IF (gl_glac(jl,jb) > 0.5_dp) THEN
!           jlat = jb
!           jlon = jl
!           xjlon = REAL(jl,dp)
!           xjlat = REAL(jb,dp)

!           !        *** search for closest Ocean Point in Ocean grid

!           CALL oclook(gl_slm, nlon, ngl, jlon, jlat, xjlon, xjlat, kfound)

!           IF (kfound == 1) THEN
!             gl_zcalv(jlon,jlat) = gl_zcalv(jlon,jlat)+gl_apmecal(jl,jb)*gridarea(jb)
!             kfocim(jlon,jlat) = 1
!             IF (ique /= 0) THEN
!               WRITE(message_text,*) '--> jlon = ', jlon, '  jlat = ', jlat
!               CALL message('glacier_to_ocean', message_text)
!             END IF
!           ELSE
!             IF (ique /= 0) THEN
!               WRITE(message_text,*) 'No ocean point at (jl,jb)=', jl,jb
!               CALL message('glacier_to_ocean', message_text)
!             END IF
!             xidb = xidb + gl_apmecal(jl,jb)
!           ENDIF
!         ENDIF
!       ENDDO
!     ENDDO

!     !    Redistributing the Water in XIDB onto all Ocean Inflow Points
!     !    Applying a weight to treat arid and humid regions differently

!     kfound = SUM(kfocim(:,:))

!     IF (kfound > 0) THEN 
!       gl_zcalv(:,:) = gl_zcalv(:,:)+gl_zcalv(:,:)/SUM(gl_zcalv(:,:)) &
!            *xidb*REAL(kfocim(:,:),dp)
!       gl_zcalv(:,:) = MAX(gl_zcalv(:,:),0.0_dp)
!     ELSE
!       WRITE(message_text,*) &
!            'Error no inflow points on ocean grid found'
!       CALL message('glacier_to_ocean', message_text)
!     ENDIF

!   END SUBROUTINE glacier_to_ocean

  SUBROUTINE glacier_to_ocean(lhd_que)
    !   ***** Glacier calving subroutine fuer ECHAM5
    !         Distributes the P-E differences over glacier grid boxes
    !         to the nearest ocean grid box indicated by the land sea mask slm.
    !         For accurate distribution, the gridbox area of source and target
    !         grid box have to be taken into account.
    !
    !    Programmed and developed by Stefan Hagemann, MPI
    !
    !   ***** Version 1.0 -- February 2001
    !         Programmierung analogous to hd_tooc.f90
    !         Subroutine uses oclook.f that is included in hd_tooc.f
    !
    !   ****** Variablenliste
    !
    !     apmecal = P minus E array -- expected in [m/s]
    !      xidb = Summation array of Inflows, for which no Inflowbox in
    !             Ocean grid was found.
    !    kfound = Inflow-Point found  on ocean grid YES/NO
    !
    !      nlon = Longitudes of global atmosphere/ocean grid
    !       ngl = Latitudes of global ocean grid
    !
    !       slm = Land Sea Mask on atmosphere/ocean grid
    USE messy_main_grid_def_mem_bi,  ONLY: nlon, ngl
    USE messy_main_grid_def_bi,      ONLY:gridarea

    LOGICAL, INTENT(in) :: lhd_que

    INTEGER :: jl, jb
    TYPE(cart_idx_2d) :: dest
    LOGICAL :: any_ocinflow
    ! mz_bk_20120802+
    ! REAL(dp) :: zcalv(nlon, ngl)
    ! mz_bk_20120802-
    REAL(dp) :: xidb

    ! mz_bk_20120802+
    ! zcalv(:,:)  = 0.0_dp
    gl_zcalv(:,:)  = 0.0_dp
    ! mz_bk_20120802-
    xidb        = 0.0_dp
    any_ocinflow = .FALSE.

    !    Redrisbution of P-E values over glacier points to closest
    !    ocean points
    !
    !   ******* Loop over all glacier points

    DO jb= 1, ngl
      DO jl= 1, nlon
        IF (gl_glac(jl,jb) > 0.5_dp) THEN
          !        *** put at closest Ocean Point in Ocean grid
          dest = oclook(gl_slm, jl, jb, cart_coord_2d(REAL(jl, dp), REAL(jb, dp)))
          IF (dest%ilon /= -1) THEN
             ! mz_bk_20120802+
            ! zcalv(dest%ilon,dest%ilat) = zcalv(dest%ilon,dest%ilat) &
            !      + gl_apmecal(jl,jb)*gridarea(jb)
            gl_zcalv(dest%ilon,dest%ilat) = gl_zcalv(dest%ilon,dest%ilat) &
                 + gl_apmecal(jl,jb)*gridarea(jb)
             ! mz_bk_20120802-
            any_ocinflow = .TRUE.
            IF (lhd_que) THEN
              WRITE(message_text,*) '--> jlon = ', dest%ilon, &
                   '  jlat = ', dest%ilat
              CALL message('glacier_to_ocean', message_text)
            END IF
          ELSE
            IF (lhd_que) THEN
              WRITE(message_text,*) 'No ocean point at (jl,jb)=', jl,jb
              CALL message('glacier_to_ocean', message_text)
            END IF
            xidb = xidb + gl_apmecal(jl,jb)*gridarea(jb)
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    !    Redistributing the Water in XIDB onto all Ocean Inflow Points
    !    Applying a weight to treat arid and humid regions differently
    ! PARALLEL: insert OR reduce on xidb /= 0.0_dp here, then sum if required
    IF (xidb /= 0.0_dp) THEN
      IF (any_ocinflow) THEN
         ! mz_bk_20120802+
        ! zcalv(:,:) = zcalv(:,:)+zcalv(:,:)/SUM(zcalv(:,:)) * xidb
        gl_zcalv(:,:) = gl_zcalv(:,:)+gl_zcalv(:,:)/SUM(gl_zcalv(:,:)) * xidb
         ! mz_bk_20120802-
      ELSE
        WRITE(message_text,*) &
             'Error no inflow points on ocean grid found'
        CALL message('glacier_to_ocean', message_text)
      ENDIF

      ! mz_bk_20120730+
      !done outside...
    ! !   Add calved water to fresh water flux into the ocean
    ! !   and subtract latent heats from heat flux and conductive heat

    !   IF (lcouple) THEN
    !     DO jb = 1, ngl
    !       DO jl = 1, nlon
    !         gl_awfre(jl,jb) = gl_awfre(jl,jb)+zcalv(jl,jb)/gridarea(jb) &
    !              /MAX(1.0_dp-gl_slf(jl,jb)-gl_alake(jl,jb),1.e-6_dp)
    !         gl_awhea(jl,jb) = gl_awhea(jl,jb)-zcalv(jl,jb)/gridarea(jb) &
    !              /MAX(1.0_dp-gl_slf(jl,jb)-gl_alake(jl,jb),1.e-6_dp)*rhoh2o*alf
    !         gl_aicon(jl,jb) = gl_aicon(jl,jb)-zcalv(jl,jb)/gridarea(jb) &
    !              /MAX(1.0_dp-gl_slf(jl,jb)-gl_alake(jl,jb),1.e-6_dp)*rhoh2o*alf
    !       ENDDO
    !     ENDDO
    !   END IF

    !   !      add 'zcalv' to 'disch' for diagnostics on 'runtoc'
    !   !      (assuming hd_model was called before !)

    !   DO jb = 1, ngl
    !     DO jl = 1, nlon
    !       gl_disch(jl,jb) = gl_disch(jl,jb)+zcalv(jl,jb)/gridarea(jb)
      ! mz_bk_20120730-
    !     ENDDO
    !   ENDDO
    END IF

    ! mz_bk_20120730+
    ! CALL glacier_diags
    ! mz_bk_20120730-

  END SUBROUTINE glacier_to_ocean
  ! mz_bk_20120730-

!-------------------------------------------------------------------------------

  ! mz_bk_20120730+
  ! SUBROUTINE oclook(foclsm, nlon, nlat, jlon, jlat, xjlon, xjlat, kfound)

  !   !*************************************************************************
  !   !
  !   ! **** Routine that looks for the next ocean gridbox closest to the
  !   !         Index jlon,jlat in the Land sea mask array foclsm(nlon,nlat)
  !   !         that corresponds to the rational index xjlon, xjlat
  !   !
  !   ! kfound = Inflow-Point found on ocean grid YES/NO
  !   !    ndd = If no Inflow-Point is found as direct neighbour,
  !   !          it is searched in NWSE direction until the maximum distance of
  !   !          NDD Boxes is reached.
  !   !          Currently it is:  ndd=INT(nlat/12): T42: 5 --> ca. 1400 km
  !   !          T106: 13 --> ca. 1430 km
  !   !          0.5 Grad: 30 --> ca. 1500 km
  !   !
  !   ! ***** Programmed and developed by Stefan Hagemann, MPI
  !   !
  !   ! ***** Version 1.0 -- Oktober 1999
  !   !
  !   INTEGER, INTENT(in)    :: nlon, nlat
  !   INTEGER, INTENT(inout) :: jlon, jlat

  !   REAL(dp) :: foclsm(nlon,nlat), xjlon, xjlat, dx, dxmin
  !   INTEGER :: ioc(4), isum, kfound

  !   INTEGER :: idd, ndd, il, ib
  !   !
  !   !  N,S,W,E-Directions
  !   ioc(:)=0
  !   IF (jlon /= 1) THEN
  !     IF (foclsm(jlon-1,jlat) < 0.5_dp) ioc(1) = 1
  !   ELSE
  !     IF (foclsm(nlon,jlat) < 0.5_dp) ioc(1) = 1
  !   ENDIF
  !   IF (jlon /= nlon) THEN
  !     IF (foclsm(jlon+1,jlat) < 0.5_dp) ioc(2) = 1
  !   ELSE
  !     IF (foclsm(1,jlat) < 0.5_dp) ioc(2) = 1
  !   ENDIF
  !   IF (jlat /= 1) THEN
  !     IF (foclsm(jlon,jlat-1) < 0.5_dp) ioc(3) = 1
  !   ENDIF
  !   IF (jlat /= nlat) THEN
  !     IF (foclsm(jlon,jlat+1) < 0.5_dp) ioc(4) = 1
  !   ENDIF
  !   isum = ioc(1)+ioc(2)+ioc(3)+ioc(4)
  !   !
  !   dxmin = 1.e9_dp
  !   IF (isum /= 0) THEN
  !     IF (ioc(1) == 1) THEN
  !       dx = SQRT( (xjlon-jlon+1)*(xjlon-jlon+1)+ &
  !                  (xjlat-jlat)*(xjlat-jlat) )
  !       IF (dx < dxmin) THEN
  !         dxmin = dx
  !         il = jlon-1
  !         ib = jlat
  !       ENDIF
  !     ENDIF
  !     IF (ioc(2) == 1) THEN
  !       dx = SQRT( (xjlon-jlon-1)*(xjlon-jlon-1)+ &
  !                  (xjlat-jlat)*(xjlat-jlat) )
  !       IF (dx < dxmin) THEN
  !         dxmin = dx
  !         il = jlon+1
  !         ib = jlat
  !       ENDIF
  !     ENDIF
  !     IF (ioc(3) == 1) THEN
  !       dx = SQRT( (xjlon-jlon)*(xjlon-jlon)+ &
  !                  (xjlat-jlat+1)*(xjlat-jlat+1) )
  !       IF (dx < dxmin) THEN
  !         dxmin = dx
  !         il = jlon
  !         ib = jlat-1
  !       ENDIF
  !     ENDIF
  !     IF (ioc(4) == 1) THEN
  !       dx = SQRT( (xjlon-jlon)*(xjlon-jlon)+ &
  !                  (xjlat-jlat-1)*(xjlat-jlat-1) )
  !       IF (dx < dxmin) THEN
  !         dxmin = dx
  !         il = jlon
  !         ib = jlat+1
  !       ENDIF
  !     ENDIF
  !     IF (il == 0) il=nlon
  !     IF (il == nlon+1) il = 1
  !     jlon = il
  !     jlat = ib
  !     kfound = 1
  !     RETURN
  !   ENDIF
  !   !
  !   !  NW,NE,SW,SE-Directions
  !   ioc(:) = 0
  !   IF (jlon /= 1) THEN
  !     IF (jlat /= 1) THEN
  !       IF (foclsm(jlon-1,jlat-1) < 0.5_dp) ioc(1) = 1
  !     ENDIF
  !     IF (jlat /= nlat) THEN
  !       IF (foclsm(jlon-1,jlat+1) < 0.5_dp) ioc(2) = 1
  !     ENDIF
  !   ELSE
  !     IF (jlat /= 1) THEN
  !       IF (foclsm(nlon,jlat-1) < 0.5_dp) ioc(1) = 1
  !     ENDIF
  !     IF (jlat /= nlat) THEN
  !       IF (foclsm(nlon,jlat+1) < 0.5_dp) ioc(2) = 1
  !     ENDIF
  !   ENDIF
  !   IF (jlon /= nlon) THEN
  !     IF (jlat /= 1) THEN
  !       IF (foclsm(jlon+1,jlat-1) < 0.5_dp) ioc(3) = 1
  !     ENDIF
  !     IF (jlat /= nlat) THEN
  !       IF (foclsm(jlon+1,jlat+1) < 0.5_dp) ioc(4) = 1
  !     ENDIF
  !   ELSE
  !     IF (jlat /= 1) THEN
  !       IF (foclsm(1,jlat-1) < 0.5_dp) ioc(3) = 1
  !     ENDIF
  !     IF (jlat /= nlat) THEN
  !       IF (foclsm(1,jlat+1) < 0.5_dp) ioc(4) = 1
  !     ENDIF
  !   ENDIF
  !   isum = ioc(1)+ioc(2)+ioc(3)+ioc(4)
  !   !
  !   dxmin = 1.e9_dp
  !   IF (isum /= 0) THEN
  !     IF (ioc(1) == 1) THEN
  !       dx = SQRT( (xjlon-jlon+1)*(xjlon-jlon+1)+ &
  !                  (xjlat-jlat+1)*(xjlat-jlat+1) )
  !       IF (dx < dxmin) THEN
  !         dxmin = dx
  !         il = jlon-1
  !         ib = jlat-1
  !       ENDIF
  !     ENDIF
  !     IF (ioc(2) == 1) THEN
  !       dx = SQRT( (xjlon-jlon+1)*(xjlon-jlon+1)+ &
  !                  (xjlat-jlat-1)*(xjlat-jlat-1) )
  !       IF (dx < dxmin) THEN
  !         dxmin = dx
  !         il = jlon-1
  !         ib = jlat+1
  !       ENDIF
  !     ENDIF
  !     IF (ioc(3) == 1) THEN
  !       dx = SQRT( (xjlon-jlon-1)*(xjlon-jlon-1)+ &
  !            &               (xjlat-jlat+1)*(xjlat-jlat+1) )
  !       IF (dx < dxmin) THEN
  !         dxmin = dx
  !         il = jlon+1
  !         ib = jlat-1
  !       ENDIF
  !     ENDIF
  !     IF (ioc(4) == 1) THEN
  !       dx = SQRT( (xjlon-jlon-1)*(xjlon-jlon-1)+ &
  !                  (xjlat-jlat-1)*(xjlat-jlat-1) )
  !       IF (dx < dxmin) THEN
  !         dxmin = dx
  !         il = jlon+1
  !         ib = jlat+1
  !       ENDIF
  !     ENDIF
  !     IF (il == 0)      il = nlon
  !     IF (il == nlon+1) il = 1
  !     jlon = il
  !     jlat = ib
  !     kfound = 1
  !     RETURN
  !   ENDIF
  !   !
  !   ! **** Second to fifth next Gridboxes
  !   !
  !   !  N,S,W,E-Directions
  !   ndd = INT(nlat/12)
  !   DO idd = 2, ndd
  !     ioc(:) = 0
  !     !
  !     !   *** Distance to central gridbox: IL, IB
  !     IF (jlon-idd >= 1) THEN
  !       IF (foclsm(jlon-idd,jlat) < 0.5_dp) ioc(1) = 1
  !     ELSE
  !       IF (foclsm(nlon+jlon-idd,jlat) < 0.5_dp) ioc(1) = 1
  !     ENDIF
  !     IF (jlon+idd <= nlon) THEN
  !       IF (foclsm(jlon+idd,jlat) < 0.5_dp) ioc(2) = 1
  !     ELSE
  !       IF (foclsm(jlon+idd-nlon,jlat) < 0.5_dp) ioc(2) = 1
  !     ENDIF
  !     IF (jlat-idd >= 1) THEN
  !       IF (foclsm(jlon,jlat-idd) < 0.5_dp) ioc(3) = 1
  !     ENDIF
  !     IF (jlat+idd <= nlat) THEN
  !       IF (foclsm(jlon,jlat+1) < 0.5_dp) ioc(4) = 1
  !     ENDIF
  !     isum = ioc(1)+ioc(2)+ioc(3)+ioc(4)
  !     !
  !     dxmin = 1.e9_dp
  !     IF (isum /= 0) THEN
  !       IF (ioc(1) == 1) THEN
  !         dx = SQRT( (xjlon-jlon+idd)*(xjlon-jlon+idd)+ &
  !                    (xjlat-jlat)*(xjlat-jlat) )
  !         IF (dx < dxmin) THEN
  !           dxmin = dx
  !           il = jlon-idd
  !           ib = jlat
  !         ENDIF
  !       ENDIF
  !       IF (ioc(2) == 1) THEN
  !         dx = SQRT( (xjlon-jlon-idd)*(xjlon-jlon-idd)+ &
  !                    (xjlat-jlat)*(xjlat-jlat) )
  !         IF (dx < dxmin) THEN
  !           dxmin = dx
  !           il = jlon+idd
  !           ib = jlat
  !         ENDIF
  !       ENDIF
  !       IF (ioc(3) == 1) THEN
  !         dx = SQRT( (xjlon-jlon)*(xjlon-jlon)+         &
  !                    (xjlat-jlat+idd)*(xjlat-jlat+idd) )
  !         IF (dx < dxmin) THEN
  !           dxmin = dx
  !           il = jlon
  !           ib = jlat-idd
  !         ENDIF
  !       ENDIF
  !       IF (ioc(4) == 1) THEN
  !         dx = SQRT( (xjlon-jlon)*(xjlon-jlon)+      &
  !                    (xjlat-jlat-idd)*(xjlat-jlat-idd) )
  !         IF (dx < dxmin) THEN
  !           dxmin = dx
  !           il = jlon
  !           ib = jlat+idd
  !         ENDIF
  !       ENDIF
  !       IF (il < 1)    il = il+nlon
  !       IF (il > nlon) il = il-nlon
  !       jlon = il
  !       jlat = ib
  !       kfound = 1
  !       RETURN
  !     ENDIF
  !     !
  !   ENDDO
  !   !
  !   kfound = 0

  ! END SUBROUTINE oclook

  PURE FUNCTION oclook(foclsm, jlon, jlat, coord) RESULT(idx)

    !*************************************************************************
    !
    ! **** Routine that looks for the next ocean gridbox closest to the
    !         Index jlon,jlat in the Land sea mask array foclsm(:,:)
    !         that corresponds to the rational index coord%lon, coord%lat
    !
    ! lfound = Inflow-Point found on ocean grid YES/NO
    !    ndd = If no Inflow-Point is found as direct neighbour,
    !          it is searched in NWSE direction until the maximum distance of
    !          NDD Boxes is reached.
    !          Currently it is:  ndd=INT(nlat/12): T42: 5 --> ca. 1400 km
    !          T106: 13 --> ca. 1430 km
    !          0.5 Grad: 30 --> ca. 1500 km
    !
    ! ***** Programmed and developed by Stefan Hagemann, MPI
    !
    ! ***** Version 1.0 -- Oktober 1999
    !
    REAL(dp), INTENT(in) :: foclsm(:,:)
    INTEGER, INTENT(in)  :: jlon, jlat
    TYPE(cart_coord_2d), INTENT(in) :: coord
    TYPE(cart_idx_2d) :: idx
    real(dp) :: dx, dxmin
    INTEGER :: ioc(4), isum

    INTEGER :: idd, ndd, il, ib, nlon, nlat
    !
    nlon = SIZE(foclsm, 1)
    nlat = SIZE(foclsm, 2)
    !
    !  N,S,W,E-Directions
    ioc(:)=0
    IF (jlon /= 1) THEN
      IF (foclsm(jlon-1,jlat) < 0.5_dp) ioc(1) = 1
    ELSE
      IF (foclsm(nlon,jlat) < 0.5_dp) ioc(1) = 1
    ENDIF
    IF (jlon /= nlon) THEN
      IF (foclsm(jlon+1,jlat) < 0.5_dp) ioc(2) = 1
    ELSE
      IF (foclsm(1,jlat) < 0.5_dp) ioc(2) = 1
    ENDIF
    IF (jlat /= 1) THEN
      IF (foclsm(jlon,jlat-1) < 0.5_dp) ioc(3) = 1
    ENDIF
    IF (jlat /= nlat) THEN
      IF (foclsm(jlon,jlat+1) < 0.5_dp) ioc(4) = 1
    ENDIF
    isum = ioc(1)+ioc(2)+ioc(3)+ioc(4)
    !
    dxmin = 1.e9_dp
    IF (isum /= 0) THEN
      IF (ioc(1) == 1) THEN
        dx = SQRT( (coord%lon-jlon+1)*(coord%lon-jlon+1)+ &
                   (coord%lat-jlat)*(coord%lat-jlat) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon-1
          ib = jlat
        ENDIF
      ENDIF
      IF (ioc(2) == 1) THEN
        dx = SQRT( (coord%lon-jlon-1)*(coord%lon-jlon-1)+ &
                   (coord%lat-jlat)*(coord%lat-jlat) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon+1
          ib = jlat
        ENDIF
      ENDIF
      IF (ioc(3) == 1) THEN
        dx = SQRT( (coord%lon-jlon)*(coord%lon-jlon)+ &
                   (coord%lat-jlat+1)*(coord%lat-jlat+1) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon
          ib = jlat-1
        ENDIF
      ENDIF
      IF (ioc(4) == 1) THEN
        dx = SQRT( (coord%lon-jlon)*(coord%lon-jlon)+ &
                   (coord%lat-jlat-1)*(coord%lat-jlat-1) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon
          ib = jlat+1
        ENDIF
      ENDIF
      IF (il == 0) il=nlon
      IF (il == nlon+1) il = 1
      idx%ilon = il
      idx%ilat = ib
      RETURN
    ENDIF
    !
    !  NW,NE,SW,SE-Directions
    ioc(:) = 0
    IF (jlon /= 1) THEN
      IF (jlat /= 1) THEN
        IF (foclsm(jlon-1,jlat-1) < 0.5_dp) ioc(1) = 1
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(jlon-1,jlat+1) < 0.5_dp) ioc(2) = 1
      ENDIF
    ELSE
      IF (jlat /= 1) THEN
        IF (foclsm(nlon,jlat-1) < 0.5_dp) ioc(1) = 1
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(nlon,jlat+1) < 0.5_dp) ioc(2) = 1
      ENDIF
    ENDIF
    IF (jlon /= nlon) THEN
      IF (jlat /= 1) THEN
        IF (foclsm(jlon+1,jlat-1) < 0.5_dp) ioc(3) = 1
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(jlon+1,jlat+1) < 0.5_dp) ioc(4) = 1
      ENDIF
    ELSE
      IF (jlat /= 1) THEN
        IF (foclsm(1,jlat-1) < 0.5_dp) ioc(3) = 1
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(1,jlat+1) < 0.5_dp) ioc(4) = 1
      ENDIF
    ENDIF
    isum = ioc(1)+ioc(2)+ioc(3)+ioc(4)
    !
    dxmin = 1.e9_dp
    IF (isum /= 0) THEN
      IF (ioc(1) == 1) THEN
        dx = SQRT( (coord%lon-jlon+1)*(coord%lon-jlon+1)+ &
                   (coord%lat-jlat+1)*(coord%lat-jlat+1) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon-1
          ib = jlat-1
        ENDIF
      ENDIF
      IF (ioc(2) == 1) THEN
        dx = SQRT( (coord%lon-jlon+1)*(coord%lon-jlon+1)+ &
                   (coord%lat-jlat-1)*(coord%lat-jlat-1) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon-1
          ib = jlat+1
        ENDIF
      ENDIF
      IF (ioc(3) == 1) THEN
        dx = SQRT( (coord%lon-jlon-1)*(coord%lon-jlon-1)+ &
             &               (coord%lat-jlat+1)*(coord%lat-jlat+1) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon+1
          ib = jlat-1
        ENDIF
      ENDIF
      IF (ioc(4) == 1) THEN
        dx = SQRT( (coord%lon-jlon-1)*(coord%lon-jlon-1)+ &
                   (coord%lat-jlat-1)*(coord%lat-jlat-1) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon+1
          ib = jlat+1
        ENDIF
      ENDIF
      IF (il == 0)      il = nlon
      IF (il == nlon+1) il = 1
      idx%ilon = il
      idx%ilat = ib
      RETURN
    ENDIF
    !
    ! **** Second to fifth next Gridboxes
    !
    !  N,S,W,E-Directions
    ndd = INT(nlat/12)
    DO idd = 2, ndd
      ioc(:) = 0
      !
      !   *** Distance to central gridbox: IL, IB
      IF (jlon-idd >= 1) THEN
        IF (foclsm(jlon-idd,jlat) < 0.5_dp) ioc(1) = 1
      ELSE
        IF (foclsm(nlon+jlon-idd,jlat) < 0.5_dp) ioc(1) = 1
      ENDIF
      IF (jlon+idd <= nlon) THEN
        IF (foclsm(jlon+idd,jlat) < 0.5_dp) ioc(2) = 1
      ELSE
        IF (foclsm(jlon+idd-nlon,jlat) < 0.5_dp) ioc(2) = 1
      ENDIF
      IF (jlat-idd >= 1) THEN
        IF (foclsm(jlon,jlat-idd) < 0.5_dp) ioc(3) = 1
      ENDIF
      IF (jlat+idd <= nlat) THEN
        IF (foclsm(jlon,jlat+1) < 0.5_dp) ioc(4) = 1
      ENDIF
      isum = ioc(1)+ioc(2)+ioc(3)+ioc(4)
      !
      dxmin = 1.e9_dp
      IF (isum /= 0) THEN
        IF (ioc(1) == 1) THEN
          dx = SQRT( (coord%lon-jlon+idd)*(coord%lon-jlon+idd)+ &
                     (coord%lat-jlat)*(coord%lat-jlat) )
          IF (dx < dxmin) THEN
            dxmin = dx
            il = jlon-idd
            ib = jlat
          ENDIF
        ENDIF
        IF (ioc(2) == 1) THEN
          dx = SQRT( (coord%lon-jlon-idd)*(coord%lon-jlon-idd)+ &
                     (coord%lat-jlat)*(coord%lat-jlat) )
          IF (dx < dxmin) THEN
            dxmin = dx
            il = jlon+idd
            ib = jlat
          ENDIF
        ENDIF
        IF (ioc(3) == 1) THEN
          dx = SQRT( (coord%lon-jlon)*(coord%lon-jlon)+         &
                     (coord%lat-jlat+idd)*(coord%lat-jlat+idd) )
          IF (dx < dxmin) THEN
            dxmin = dx
            il = jlon
            ib = jlat-idd
          ENDIF
        ENDIF
        IF (ioc(4) == 1) THEN
          dx = SQRT( (coord%lon-jlon)*(coord%lon - jlon)+      &
                     (coord%lat-jlat-idd)*(coord%lat-jlat-idd) )
          IF (dx < dxmin) THEN
            dxmin = dx
            il = jlon
            ib = jlat+idd
          ENDIF
        ENDIF
        IF (il < 1)    il = il+nlon
        IF (il > nlon) il = il-nlon
        idx%ilon = il
        idx%ilat = ib
        RETURN
      ENDIF
      !
    ENDDO
    !
    idx%ilon = -1
    idx%ilat = -1
  END FUNCTION oclook

  SUBROUTINE hydrology_slm_invariants
    USE messy_main_mpi_bi,             ONLY: gather_field
    USE messy_main_grid_def_mem_bi,    ONLY: nlon, ngl
    USE messy_main_data_bi,            ONLY: slm, slf, alake
    REAL(dp) :: joined_slm(nlon, ngl)
    REAL(dp), POINTER :: gl(:,:)
    gl => gl_slm
!    CALL gather_gp (gl, slm, dcg)
    CALL gather_field(gl, slm)
    gl => gl_slf
!    CALL gather_gp (gl, slf, dcg)
    CALL gather_field(gl, slf)
    gl => gl_alake
!    CALL gather_gp (gl, alake, dcg)
    CALL gather_field(gl, alake)

    CALL slm_lake_join(ngl, nlon, joined_slm, gl_slm, gl_slf, gl_alake)

    IF (p_pe == p_io) THEN
      CALL fill_oclook_caches(joined_slm)
      CALL prepare_intpol_mapping
      CALL create_kasglob_list(alf_n, alf_k, overlandflow, alf_n_kas)
      CALL create_kasglob_list(arf_n, arf_k, riverflow, arf_n_kas)
    END IF
  END SUBROUTINE hydrology_slm_invariants

  ! mz_bk_20120730+
  ! SUBROUTINE cleanup_hydrology_slm_invariants
  ! mz_bk_20120730-
  SUBROUTINE cleanup_hydrology_slm_invar
    IF (p_pe == p_io) THEN
      CALL cleanup_kasglob_list(alf_n_kas)
      CALL cleanup_kasglob_list(arf_n_kas)
    END IF
    ! mz_bk_20120730+
  ! END SUBROUTINE cleanup_hydrology_slm_invariants
    ! mz_bk_20120730-
  END SUBROUTINE cleanup_hydrology_slm_invar

  SUBROUTINE fill_oclook_caches(slm)
    USE messy_main_grid_def_mem_bi,  ONLY: nlon, ngl
    USE messy_main_grid_def_bi,      ONLY:philat
    REAL(dp), INTENT(in) :: slm(nlon, ngl)

    REAL(dp) :: fb, fl, xjlon, xanf, xend, ocscal
    INTEGER :: jlon, jlat, jb, jl
    TYPE(cart_coord_2d) :: coord

    ocscal = fullcirc / REAL(nlon, dp)
    DO jb = 1, nb
      fb = fborg-REAL(jb,dp)*fscal+0.5_dp*fscal
      ! mz_bk_20120730+
      ! jlat = dec_monotonic_interval_closest_midpoint(philat, &
      !      fb, aub=90._dp, alb=-90._dp)
      ! mz_bk_20120730-
      jlat = dec_mon_interval_clos_midpoint(philat, &
           fb, aub=90._dp, alb=-90._dp)
      IF (jlat == 1) THEN
        xanf = 90.0_dp
        xend = (philat(2)+philat(1))*0.5_dp
      ELSEIF (jlat == ngl) THEN
        xanf = (philat(ngl)+philat(ngl - 1))*0.5_dp
        xend = -90.0_dp
      ELSE
        xanf = (philat(jlat - 1) + philat(jlat))*0.5_dp
        xend = (philat(jlat + 1) + philat(jlat))*0.5_dp
      ENDIF
      coord%lat = REAL(jlat, dp) + (xanf-fb) / (xanf-xend)
      DO jl = 1, nl
        fl = MOD(REAL(jl,dp)*fscal+florg-0.5_dp*fscal + 180._dp, &
             fullcirc) - 180.0_dp
        xjlon = (fl - oclorg + ocscal * 0.5_dp)/ocscal + 1.0_dp
        IF (INT(xjlon + 0.00001_dp) <= 0) xjlon = xjlon + REAL(nlon, dp)
        jlon = INT(xjlon + 0.00001_dp)
        coord%lon = xjlon
        oclook_cache(jl, jb) = oclook(slm, jlon, jlat, coord)
        IF (oclook_cache(jl, jb)%ilat == -1 .AND. lhd_que) THEN
          WRITE(message_text,*) 'potential problem: no inflow points on &
               & ocean grid found for idx(', jl, ',', jb, ')'
          CALL message('hd: fill_oclook_caches', message_text)
        END IF
      END DO
    END DO
  END SUBROUTINE fill_oclook_caches

  SUBROUTINE create_kasglob_list(a_n, a_k, iflow, a_n_kas)
    REAL(dp), INTENT(in) :: a_n(:, :), a_k(:, :)
    INTEGER, INTENT(in) :: iflow
    TYPE(cart_xidx_2d), ALLOCATABLE, INTENT(inout) :: a_n_kas(:)
    REAL(dp) :: divmm
    INTEGER :: size_i, size_j, i, j, num_extents, extent, extlen
    INTEGER :: jl, extelem

    IF (iflow == overlandflow) THEN
      divmm = 1.0_dp
    ELSE IF (iflow == riverflow) THEN
      divmm = 1.0_dp/REAL(mm,dp)
    ENDIF

    size_i = SIZE(a_n, 1)
    size_j = SIZE(a_n, 2)
    num_extents = 0
    IF (ALLOCATED(a_n_kas)) THEN
      CALL cleanup_kasglob_list(a_n_kas)
    END IF
    DO j = 1, size_j
      i = 1
      DO WHILE(i <= size_i)
        IF (a_n(i, j) > 0.5_dp) THEN
          num_extents = num_extents + 1
          i = i + 1
          DO WHILE(i <= size_i)
            IF (a_n(i, j) <= 0.5_dp) EXIT
            i = i + 1
          END DO
        ELSE
          i = i + 1
        END IF
      END DO
    END DO
    ALLOCATE(a_n_kas(num_extents))
    extent = 0
    DO j = 1, size_j
      i = 1
      DO WHILE(i <= size_i)
        IF (a_n(i, j) > 0.5_dp) THEN
          extent = extent + 1
          extlen = 0
          a_n_kas(extent)%ilat = j
          a_n_kas(extent)%ilon = i
          DO WHILE(i <= size_i)
            IF (a_n(i, j) <= 0.5_dp) EXIT
            i = i + 1
            extlen = extlen + 1
          END DO
          a_n_kas(extent)%extlen = extlen
          ALLOCATE(a_n_kas(extent)%amod(extlen), a_n_kas(extent)%akdiv(extlen))
          DO extelem = 1, extlen
            jl = a_n_kas(extent)%ilon + extelem - 1
            a_n_kas(extent)%amod(extelem) = a_k(jl,j) * a_n(jl,j) &
                 / AINT(a_n(jl,j))
            a_n_kas(extent)%akdiv(extelem) = 1.0_dp &
                 / (a_n_kas(extent)%amod(extelem) + divmm)
          END DO
        ELSE
          i = i + 1
        END IF
      END DO
    END DO
  END SUBROUTINE create_kasglob_list

  SUBROUTINE cleanup_kasglob_list(a_n_kas)
    TYPE(cart_xidx_2d), ALLOCATABLE, INTENT(inout) :: a_n_kas(:)
    INTEGER :: i, n
    n = SIZE(a_n_kas)
    DO i = 1, n
      DEALLOCATE(a_n_kas(i)%amod, a_n_kas(i)%akdiv)
    END DO
    DEALLOCATE(a_n_kas)
  END SUBROUTINE cleanup_kasglob_list

  SUBROUTINE prepare_intpol_mapping
    USE messy_main_grid_def_mem_bi,            ONLY: nlon, ngl
    REAL(dp) :: axlon(nlon + 1), axlat(ngl + 2), xr(nl), yr(nb)
    CALL intpol_coord_axis_setup(axlon, axlat, xr, yr)
    CALL intpol_compute_mapping(intpol_mapping, &
         nlon + 1, ngl + 2, axlon, axlat, nl, nb, xr, yr)
  END SUBROUTINE prepare_intpol_mapping

  ! mz_bk_20120730-

  ! mz_bk_20120730+
  ! from mo_array_utils.f90 (ECHAM6/MPIESM-1.0.00)
  !> search for i such that (a(i - 1) + a(i))/2 >= x > (a(i) + a(i + 1))/2
  !> a satisfies a(i) >= a(i + 1)
  !> also considers aub >= x > (a(1) + a(2))/2
  !> and (a(size(a) - 1) + a(size(a)))/2 >= x > alb
  !> if alb or aub are given respectively
  !> returns -1 if no index satisfies the condition
  !PURE FUNCTION dec_monotonic_interval_closest_midpoint(a, x, alb, aub, first_i) RESULT(i)
  PURE FUNCTION dec_mon_interval_clos_midpoint(a, x, alb, aub, first_i) RESULT(i)
    REAL(dp), INTENT(in)           :: a(:), x
    REAL(dp), INTENT(in), OPTIONAL :: alb, aub
    INTEGER,  INTENT(in), OPTIONAL :: first_i

    INTEGER :: i, m, n

    n = SIZE(a)
    IF (n < 2) THEN
      i = -1
      RETURN
    END IF
    IF (PRESENT(aub)) THEN
      IF (aub >= x .AND. x > (a(1) + a(2)) * 0.5_dp) THEN
        i = 1
        RETURN
      ELSE IF(x > (a(1) + a(2)) * 0.5_dp) THEN
        i = -1
        RETURN
      END IF
    ELSE IF(x > (a(1) + a(2)) * 0.5_dp) THEN
      i = -1
      RETURN
    END IF
    IF (PRESENT(alb)) THEN
      IF ((a(n - 1) + a(n)) * 0.5_dp >= x .AND. x > alb) THEN
        i = n
        RETURN
      ELSE IF((a(n - 1) + a(n)) * 0.5_dp >= x) THEN
        i = -n
        RETURN
      END IF
    ELSE IF((a(n - 1) + a(n)) * 0.5_dp >= x) THEN
      i = -n
      RETURN
    END IF
    IF (n < 3) THEN
      i =  -1
      RETURN
    END IF
    ! at this point the following holds:
    ! 1. (a(1) + a(2))*0.5_dp < x and (a(n - 1) + a(n))*0.5_dp >= x
    ! 2. n > 2
    ! therefore it is guaranteed, we can find an index satisfying
    ! above condition
    m = 1
    i = (n + 1)/ 2
    IF (PRESENT(first_i)) THEN
      IF (first_i < n .AND. first_i > 1) i = first_i
    END IF
    DO
      IF ((a(i - 1) + a(i))*0.5_dp < x) THEN
        n = i
      ELSE IF (x <= (a(i) + a(i + 1))*0.5_dp) THEN
        m = i
      ELSE ! IF ((a(i - 1) + a(i))*0.5_dp >= x .AND. x > (a(i) + a(i + 1))/2)
        RETURN
      END IF
      i = (m + n + 1)/2
    END DO
!  END FUNCTION dec_monotonic_interval_closest_midpoint
  END FUNCTION dec_mon_interval_clos_midpoint
  ! mz_bk_20120730-

!-------------------------------------------------------------------------

  SUBROUTINE hd_read_nml_cpl(status, iou)

    ! a2o MODULE ROUTINE (ECHAM5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Pozzer Andrea, MPICH, Dec 2007

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'hd_read_nml_cpl'

    ! LOCAL
    LOGICAL         :: lex      ! file exists ?
    INTEGER         :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE hd_read_nml_cpl
!-----------------------------------------------------------------------

END MODULE messy_hd_e5
