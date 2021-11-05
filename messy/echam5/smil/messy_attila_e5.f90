! *********************************************************************
! SMIL INTERFACE TO ECHAM5 FOR ATTILA
! Version: see "modver" in messy_attila.f90
!
! Authors : Michael Traub, Patrick Joeckel, MPICH, June-Oct  2003
!           Patrick Joeckel, MPICH, DLR, 2003- ...
!           Sabine Brinkop,  MIM, DLR,   2007- ...
!
!
! NOTES:
!  - ATTILA_INTEGRATE MUST NOT be called within region-loop !!!
!    (main entry points must be global_start or global_end!).
!
! TODO:
!  - test ADICO(:) = (0,?,0) etc.
!  - move 'trajectory-mode' completely to SMIL-level ...?
!  - attila_tools: optmise MPI-routines (e.g.
!    gather_all instead of gather/bcast)
!  - clean up comments, most of them are outdated!
!  - clean up output to log-file
!  - allow chunk-size with rest
!  - better random number handling for convection
!  - test with MESSYTENDENCY
!  - check new vertical grids / velocities with and wihtout convection
!  - ...
! **********************************************************************

! **********************************************************************
! NOTES FOT THE TRAJECTORY-mode:
!     - attila_initialize_positions_t
!        1. All trajectories get their individual start-date in the namelist
!           /TRAJ/.
!        2. If LTRAJEC_SAME_DATE is T, the individual dates are overwritten
!           with LTRAJEC_DATE (during initialization !).
!       ADVANTAGE: Only one init-positions-algorithm is needed in the
!                  code (namley the one for individual trajectories.
!        3. As long as the start-date of a trajectory is not yet reached,
!           it is reset to the initial position in every time step.
!     - The latitude in the namelist is now -90 (SP) to +90 (NP), which
!       is much more user-friendly!
!       (A longitude between -180 and 180 in the namelist would be a
!       'nice to have' ... but not as straightforward ...)
! **********************************************************************

! **********************************************************************
MODULE messy_attila_e5
! **********************************************************************

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,     ONLY: start_message_bi, end_message_bi &
                                     , error_bi
  USE messy_main_rnd_bi,         ONLY: RND_F90, RND_MTW, RND_LUX      &
                                     , RND_MP_SEQ, RND_MP_PIN         &
                                     , RND_MP_PSJ, RND_MP_PIJ         &
                                     , rnd_init_bi, rnd_finish_bi     &
                                     , rnd_number_bi

  ! MESSy
  USE messy_main_constants_mem,  ONLY: dp
  USE messy_main_tools,          ONLY: PTR_2D_ARRAY
  USE messy_main_channel,        ONLY: STRLEN_OBJECT
  USE messy_attila_mem_e5                                ! op_pj_20140212
  USE messy_attila

  IMPLICIT NONE
  PRIVATE

  ! INTRINSIC PROCEDURES
  INTRINSIC :: ASSOCIATED, SIZE, EXP, INT, REAL, TRIM &
             , ALL, NULL, SUM, NINT, MAX, MIN, LOG

  ! POINTERS TO NEW CHANNEL OBJECTS
  ! LATITUDES (0 (NP)... 180 (SP)) OF CELLS
  REAL(dp), DIMENSION(:), POINTER :: latz_math => NULL()
  ! LATITUDES (-90 (SP) ... 90 (NP)) OF CELLS
  REAL(dp), DIMENSION(:), POINTER :: latz      => NULL()
  ! LONGITUDES (0 ... 360) OF CELLS
  REAL(dp), DIMENSION(:), POINTER :: lonz      => NULL()
  REAL(dp), DIMENSION(:), POINTER :: etaheight => NULL() ! ETA_HEIGHT
  REAL(dp), DIMENSION(:), POINTER :: pressheight=>NULL()! PRESSURE HEIGHT IN Pa
  REAL(dp), DIMENSION(:), POINTER :: press_m1 => NULL() ! PRESSURE HEIGHT AT TM1

  ! RANDOM NUMBER HANDLING --------------------------------------------------
  TYPE PTR_2D_ARRAY_I
     INTEGER, DIMENSION(:,:), POINTER  :: PTR => NULL()
  END TYPE PTR_2D_ARRAY_I
  !
  ! INDEX OF RANDOM NUMBER SET FOR DIFFERENT PROCESSES
  INTEGER, PARAMETER :: RANDI_TURB = 1 ! LENGTH: NCELL x 1
  INTEGER, PARAMETER :: RANDI_CONV = 2 ! LENGTH: NCELL x 1
  INTEGER, PARAMETER :: RANDI_CAT  = 3 ! LENGTH: NCELL x 1
  INTEGER, PARAMETER :: RANDI_MOCA = 4 ! LENGTH: NCELL x 3
  INTEGER, PARAMETER :: NRPROC     = 4 ! NUMBER OF RANDOM PROCESSES
  !
  ! RAND_MULT x NCELL RANDOM NUMBERS FOR EACH RANDOM PROCESS
  INTEGER, PARAMETER, DIMENSION(NRPROC) :: RAND_MULT = (/ 1, 1, 1, 3 /)
  !
  ! SWITCH FOR PROCESS
  LOGICAL,            DIMENSION(NRPROC) :: LPROC = &
       (/ .FALSE., .FALSE., .FALSE., .FALSE. /)
  !
  ! DEFAULT RANDOM NUMBER (FOR PROCESS OFF)
  REAL(DP), PARAMETER, DIMENSION(NRPROC) :: RDEF = &
       (/ 0._dp, 0._dp, 0._dp, 0.5_dp /)
  !
  ! START SEEDS
  INTEGER, PARAMETER, DIMENSION(NRPROC) :: START_SEED = &
       (/211169, 7249, 13372, 159941/)
  !
  ! WORKSPACE FOR RANDOM NUMBERS
  TYPE(PTR_2D_ARRAY),   DIMENSION(NRPROC), SAVE :: HARVEST
  !
  ! PSEUDO RANDOM NUMBER STREAM IDs
  INTEGER, DIMENSION(NRPROC), SAVE :: rndid = (/ 0, 0, 0, 0 /)
  ! PSEUDO RANDOM NUMBER GENERATOR
  INTEGER, SAVE              :: RND_ATTILA = RND_MTW  ! default
  !
  INTEGER, PARAMETER, PUBLIC :: I_RANDOM_F90 = 0  ! F90 INTRINSIC
  INTEGER, PARAMETER, PUBLIC :: I_RANDOM_MTW = 1  ! MERSENNE TWISTER
  INTEGER, PARAMETER, PUBLIC :: I_RANDOM_LUX = 2  ! LUXURY
  !
  REAL(DP), DIMENSION(:), POINTER :: test_rnd => NULL()
  ! ---------------------------------------------------------------------

  ! INDEX BOUNDARIES OF NCELL IN NGCELL
  INTEGER, DIMENSION(:,:), POINTER, SAVE :: IDX => NULL()
  INTEGER, SAVE :: np
  INTEGER, DIMENSION(:),   POINTER, SAVE :: V_NCELL => NULL()

  ! DIAGNOSTIC FIELDS IN ATTILA CHANNEL
! op_pj_20140212+
! moved to messy_attila_mem_e5 to resolve circular dependency
!!$  ! NUMBER OF CELLS PER GRID BOX (ALL CPUs) IN DECOMPOSITION
!!$  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER, SAVE :: SGNCB   => NULL()
!!$  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER, SAVE :: SGNCBM1 => NULL()
! op_pj_20140212-
  ! NUMBER OF CELLS PER GRID BOX (ALL CPUs) IN DECOMPOSITION
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER, SAVE :: SGNCBL => NULL()

  ! LG-CHANNEL OBJECTS FOR CONVECTION
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER, SAVE :: SGNCB_CUD    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER, SAVE :: SGNCB_CSUB   => NULL()
  ! ... level index change by convection
  REAL(dp),     DIMENSION(:),     POINTER, SAVE :: CONVDELK => NULL()
  ! ... number of iterations for convection
  REAL(dp),                       POINTER, SAVE :: NCONV_ITER => NULL()
  ! ... level indices before and after updraft, downdraft, subsidence
  REAL(dp), DIMENSION(:),     POINTER, SAVE :: pos_beg_up   => NULL()
  REAL(dp), DIMENSION(:),     POINTER, SAVE :: pos_end_up   => NULL()
  REAL(dp), DIMENSION(:),     POINTER, SAVE :: pos_end_do   => NULL()
  REAL(dp), DIMENSION(:),     POINTER, SAVE :: pos_end_sub  => NULL()
  REAL(DP),                   POINTER, SAVE :: rlg_conv     => NULL()

  ! ADDITIONAL VARIABLES USED WITHIN TIME LOOP (RESET AT FIRST STEP)
  LOGICAL, SAVE :: LLCONV_SAVE
  LOGICAL, SAVE :: LLCAT_SAVE
! op_sb_20140827+
  ! .TRUE., if initialisation has been finished, this happens if
  ! lstart directly before the call to integrate in global_end
  LOGICAL, SAVE :: L_INI_FINI = .FALSE.
! op_sb_20140827-
  ! LOCAL FIELDS
  REAL(dp), DIMENSION(:,:),   POINTER :: zaps_scb => NULL() ! surface pressure
  REAL(dp), DIMENSION(:,:,:), POINTER :: zpi_scb  => NULL() ! press. at interf.

  ! -> FOR ALTERNATIVE GRIDS
  ! op_sb_20140121+
  ! --> calculated here
  REAL(dp), DIMENSION(:,:,:), POINTER :: xidot        => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: xidot_cor    => NULL()
  ! diabatic temperature tendencies
#ifdef MESSYTENDENCY
  REAL(DP), DIMENSION(:,:,:), POINTER :: tte_dyn => NULL()
#endif
  REAL(DP), DIMENSION(:,:,:), POINTER :: tempdot => NULL()
  ! (theta) diabatic temperature tendencies
  REAL(DP), DIMENSION(:,:,:), POINTER :: dthetadt => NULL()
  ! --> via get_channel_object from ECHAM5 (messy_main_data_bi.f90, dyn.f90)
  REAL(dp), DIMENSION(:,:,:),   POINTER :: sigmadot  => NULL()
  REAL(dp), DIMENSION(:,:,:),   POINTER :: sidot     => NULL() ! op_sb_20160707
  REAL(dp), DIMENSION(:,:),     POINTER :: tp  => NULL()
  ! op_sb-20140602+
  REAL(dp), DIMENSION(:,:,:),   POINTER :: diab => NULL()
  REAL(dp), DIMENSION(:,:,:),   POINTER :: kine => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: sigma_ref => NULL() ! op_sb_20151113
  REAL(dp), DIMENSION(:,:),     POINTER :: i_ref => NULL()     ! op_sb_20160708

  ! op_sb-20140602-

  ! -> ONLY FOR INTERNAL PBLH CALCULATION
  REAL(dp), DIMENSION(:,:,:), POINTER :: geop_3d    => NULL() ! surf. geopot.
  REAL(dp), DIMENSION(:,:,:), POINTER :: pteta1_3d  => NULL()
  ! richardson number
  REAL(dp), DIMENSION(:,:,:), POINTER :: priad_3d   => NULL()

  ! VIA channel defined in CPL-namelist
  ! -> ONLY FOR EXTERNAL PBLH CALCULATION
  !  - planetary boundary layer height index
  REAL(dp), DIMENSION(:,:),   POINTER :: pblh_idx_gl => NULL()
  ! VIA messy_main_data_bi
  REAL(dp), DIMENSION(:,:),   POINTER :: sqcst  => NULL()    ! cos(latitude)
! op_pj_20141017+

  ! FOR CONVECTION
  REAL(dp), DIMENSION(:,:,:), POINTER :: umassfl => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: dmassfl => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: typeconv=> NULL() ! mz_mt_20041208

  ! AUXILIARY POINTERS
  ! FOR MORE INFORMATION SEE SR *ATTILA_INIT_MEMORY*
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: LSPAC      => NULL() ! APOS
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: LSPACINT   => NULL() ! NPOS
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: LSPACINTM1 => NULL() ! NPOSM1
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: LSPACPRESSM1  => NULL() ! APOSM1
  !
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: ptrlat    => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: ptrlon    => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: ptreta    => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: ptrpress  => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: ptrpressm1=> NULL()
  REAL(dp),                       POINTER :: ptramcell => NULL()
  !
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: ptrintlat => NULL()  ! ATTILA-LAT
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: ptrintlon => NULL()  ! ATTILA-LON
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: ptrintlev => NULL()  ! ATTILA-LEV
  ! ECHAM5-LON
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: ptrintlon_gp => NULL()
  ! ECHAM5-LAT
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: ptrintlat_gp => NULL()
  ! ... AT t-1
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: ptrintlatm1    => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: ptrintlonm1    => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: ptrintlevm1    => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: ptrintlonm1_gp => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: ptrintlatm1_gp => NULL()

  ! PARAMETER DIFFERENT ON EACH PROCESSOR
  ! NUMBER OF CELLS PER GRID BOX IN CORRECT DIMENSION
  REAL(dp), DIMENSION(:,:,:), POINTER :: NCB_P  => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: NCBL_P => NULL()

  ! FOR TRAJECTORY MODE ========================================
  !   TYPE FOR INITIAL POSITION I/O
  TYPE POSITION
     REAL(dp) :: lat, lon, press
     INTEGER  :: year, month, day, hour, minute
  END TYPE POSITION
  !
  INTEGER, PARAMETER                     :: MAXNAIRPARC = 300000
  INTEGER                                :: NAIRPARC
  TYPE(POSITION), DIMENSION(MAXNAIRPARC) :: AP
  ! INITIALIZATION FLAG FOR TRAJECTORY
  REAL(dp), DIMENSION(:), POINTER        :: RTRAJEC_INTEGRATE => NULL()
  ! ===========================================================

  ! CPL-namleist VARIABLES
  LOGICAL :: L_INI_PARALLEL = .TRUE.        ! initialisation in paralell
  INTEGER :: I_RANDOM_METHOD = I_RANDOM_F90 ! WHICH RANDOM NUMBER GENERATOR
  !                                         ! (default: Fortran intrinsic)
  INTEGER :: I_RANDOM_PARALLEL = RND_MP_PSJ ! HOW TO PRALLELIZE
  !                                         ! (default: sync. by jump ahead)
  LOGICAL :: L_RANDOM_TEST = .FALSE.        ! create channel object
  !                                         ! with random numbers for testing
  !
  ! channel,object for pbl-height index
  CHARACTER(LEN=STRLEN_OBJECT), PUBLIC  :: C_PBLH_INDEX(2) = ''
  ! channel,object for convection fluxes
  CHARACTER(LEN=STRLEN_OBJECT), PUBLIC  :: C_CONV_uflx(2)  = ''
  ! channel,object for convection fluxes
  CHARACTER(LEN=STRLEN_OBJECT), PUBLIC  :: C_CONV_dflx(2)  = ''
  ! channel,object for convection fluxes
  CHARACTER(LEN=STRLEN_OBJECT), PUBLIC  :: C_CONV_entr(2)  = ''
  ! channel,object for convection fluxes
  CHARACTER(LEN=STRLEN_OBJECT), PUBLIC  :: C_CONV_dntr(2)  = ''
  ! channel,object for convection fluxes
  CHARACTER(LEN=STRLEN_OBJECT), PUBLIC  :: C_CONV_tpcv(2)  = ''
  ! channel,object for convection fluxes
  CHARACTER(LEN=STRLEN_OBJECT), PUBLIC  :: C_CONV_type(2)  = ''
  ! channel,object for CAT index
  CHARACTER(LEN=STRLEN_OBJECT), PUBLIC  :: C_TURB_CATI(2)  = ''

  ! POINTERS TO 'CPL'-CHANNEL OBJECTS
  ! PBL-height index
  REAL(dp),  DIMENSION(:,:),     POINTER :: pblh_i    => NULL()
  ! upward massflux
  REAL(dp),  DIMENSION(:,:,:),   POINTER :: uflx_conv => NULL()
  ! downward massflux
  REAL(dp),  DIMENSION(:,:,:),   POINTER :: dflx_conv => NULL()
  ! entrainment
  REAL(dp),  DIMENSION(:,:,:),   POINTER :: entr_conv => NULL()
  ! detrainment
  REAL(dp),  DIMENSION(:,:,:),   POINTER :: detr_conv => NULL()
  !
  REAL(dp),  DIMENSION(:,:),     POINTER :: tpcv_conv => NULL()
  !
  REAL(dp),  DIMENSION(:,:),     POINTER :: type_conv => NULL()
  ! CAT index
  REAL(dp),  DIMENSION(:,:,:),   POINTER :: cati      => NULL()

  ! PUBLIC INTERFACE ROUTINES:
  ! ----------------------------------
  PUBLIC :: attila_initialize
  PUBLIC :: attila_init_memory
  PUBLIC :: attila_init_coupling
  PUBLIC :: attila_init_tracer
  PUBLIC :: attila_local_start
  PUBLIC :: attila_local_end
!!$  PUBLIC :: attila_global_start
  PUBLIC :: attila_global_end
  PUBLIC :: attila_free_memory

  ! PRIVATE ROUTINES:
  ! -----------------------------------
  !PRIVATE attile_read_nml_cpl
  !PRIVATE attila_initialize_positions
  !PRIVATE attila_initialize_positions_t  ! for trajectory mode
  !PRIVATE attila_integrate
  !PRIVATE attila_update_celldist
  !PRIVATE attila_convect
  !PRIVATE attila_calc_sidot
  !PRIVATE attila_calc_xidot
CONTAINS

!###########################################################################
! ### PUBLIC SUBROUTINES ###################################################
!###########################################################################

!--------------------------------------------------------------------------
  SUBROUTINE attila_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,             ONLY: p_parallel_io, p_io, p_pe, p_bcast
    USE messy_main_transform_bi,       ONLY: get_dc_index
    USE messy_main_timer,              ONLY: delta_time
    USE messy_main_grid_def_bi,        ONLY: gl_gmu
    USE messy_main_grid_def_mem_bi,    ONLY: NPLVP1, NLMSGL, NPLVP2, NLMSLP  &
                                           , apsurf, apzero, vct, nvclev, nn &
                                           , nlat=>ngl, nlon, nlev
    USE messy_main_channel_error_bi,   ONLY: channel_halt
    USE messy_main_channel_bi,         ONLY: LG_ATTILA, DC_BC, DC_IX
    USE messy_main_channel_dimensions, ONLY: new_dimension
    USE messy_main_channel_repr,       ONLY: new_representation, AUTO   &
                                           , set_representation_decomp  &
                                           , IRANK, PIOTYPE_COL
    ! MESSy
    USE messy_main_constants_mem,  ONLY: radius_earth            &
                                       , g, R_gas, M_air, cp_air
    USE messy_main_tools,          ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    REAL(DP), PARAMETER :: rd = 1000._DP * R_gas / M_air ! J/mol/kg
    CHARACTER(LEN=*), PARAMETER :: substr = 'attila_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i
    INTEGER                     :: DIMID_NGCELL
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    CALL start_message_bi(modstr, 'INITIALIZATION', substr)

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL attila_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi('error in attila_read_nml_ctrl', substr)
    END IF

    ! BROADCAST RESULTS
    ! GLOBAL SETTINGS
    CALL p_bcast(NCHUNK,          p_io)
    CALL p_bcast(CPGBAVE,         p_io)
    CALL p_bcast(LLTINFO,         p_io)
    CALL p_bcast(I_PBLH_METHOD,   p_io)
    ! PROCESS SPECIFIC SETTINGS
    CALL p_bcast(ADICO,     p_io)
    CALL p_bcast(LLTBLTURB, p_io)
    CALL p_bcast(LLCONV,    p_io)
    CALL p_bcast(LLCAT,     p_io)
    CALL p_bcast(LVDIAG,    p_io)
    ! SPECIAL MODI
    ! - RESOLUTION INDEPENDENT NUMBER OF CELLS
    CALL p_bcast(I_NCELL, p_io)
    ! - TRAJECTORY MODE
    CALL p_bcast(LTRAJEC,      p_io)
    CALL p_bcast(LTRAJEC_DATE, p_io)
    CALL p_bcast(LTRAJEC_SAME_DATE, p_io)
    ! op_sb_20140602+
    CALL p_bcast(I_VERT, p_io)
    CALL p_bcast(PRESS_REF, p_io)
    ! op_sb_20140602-
    ! SAVE VALUES
    LLCONV_SAVE = LLCONV
    LLCAT_SAVE = LLCAT

    IF (p_parallel_io) CALL ATTILA_MESSAGE

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL attila_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('error in attila_read_nml_cpl', substr)
    END IF
    !
    ! CPL NAMELIST
    CALL p_bcast(L_INI_PARALLEL, p_io)
    CALL p_bcast(I_RANDOM_METHOD, p_io)
    CALL p_bcast(I_RANDOM_PARALLEL, p_io)
    CALL p_bcast(L_RANDOM_TEST, p_io)
    ! FOR TURBULENT MIXING IN BOUNDARY LAYER
    CALL p_bcast(C_PBLH_INDEX(1), p_io)
    CALL p_bcast(C_PBLH_INDEX(2), p_io)
    ! FOR CONVECTION
    CALL p_bcast(C_CONV_uflx(1), p_io)
    CALL p_bcast(C_CONV_uflx(2), p_io)
    CALL p_bcast(C_CONV_dflx(1), p_io)
    CALL p_bcast(C_CONV_dflx(2), p_io)
    CALL p_bcast(C_CONV_entr(1), p_io)
    CALL p_bcast(C_CONV_entr(2), p_io)
    CALL p_bcast(C_CONV_dntr(1), p_io)
    CALL p_bcast(C_CONV_dntr(2), p_io)
    CALL p_bcast(C_CONV_tpcv(1), p_io)
    CALL p_bcast(C_CONV_tpcv(2), p_io)
    CALL p_bcast(C_CONV_type(1), p_io)
    CALL p_bcast(C_CONV_type(2), p_io)
    ! FOR CLEAR AIR TURBULENCE
    CALL p_bcast(C_TURB_CATI(1), p_io)
    CALL p_bcast(C_TURB_CATI(2), p_io)

    ! TRAJECTORY MODE
    IF (LTRAJEC) THEN
       CALL p_bcast(NAIRPARC, p_io)
       !
       DO i=1, NAIRPARC
          CALL p_bcast(AP(i)%lat, p_io)
          CALL p_bcast(AP(i)%lon, p_io)
          CALL p_bcast(AP(i)%press, p_io)
          CALL p_bcast(AP(i)%year, p_io)
          CALL p_bcast(AP(i)%month, p_io)
          CALL p_bcast(AP(i)%day, p_io)
          CALL p_bcast(AP(i)%hour, p_io)
          CALL p_bcast(AP(i)%minute, p_io)
       END DO
    END IF

    CALL ATTILA_GLOBAL_INIT(nlat, nlon, nlev, delta_time  &
         , .TRUE.  & ! ECHAM5 ALWAYS IN 'PARALLEL' MODE
         , vct(:), nvclev, apsurf     &
         , apzero, gl_gmu, nn, NPLVP1, NLMSGL, NPLVP2, NLMSLP, status)
    IF (status /= 0) CALL error_bi( &
         'ATTILA_GLOBAL_INIT reported an error !', substr)

    CALL ATTILA_INICOM_1

    ! SET PSEUDO RANDOM NUMBER METHOD (FOR ALL PROCESSES)
    SELECT CASE(I_RANDOM_METHOD)
    CASE(I_RANDOM_F90)
       IF (p_parallel_io) &
            WRITE(*,*) 'RANDOM NUMBERS: F90-INTRINSIC    (MACHINE DEPENDENT!)'
       RND_ATTILA = RND_F90
    CASE(I_RANDOM_MTW)
       IF (p_parallel_io) &
            WRITE(*,*) 'RANDOM NUMBERS: MERSENNE TWISTER (MACHINE INDEPENDENT!)'
       RND_ATTILA = RND_MTW
    CASE(I_RANDOM_LUX)
       IF (p_parallel_io) &
            WRITE(*,*) 'RANDOM NUMBERS: LUXURY (MACHINE INDEPENDENT!)'
       RND_ATTILA = RND_LUX
    CASE DEFAULT
       CALL error_bi('UNKNOW RANDOM NUMBER METHOD', substr)
    END SELECT

    ! LIMIT NUMBER OF CELLS FOR HIGH (MA) RESOLUTIONS (NAMELIST)
    IF (I_NCELL > 0) NGCELL = I_NCELL

    ! SET NGCELL = NAIRPARC IN TRAJECTORY MODE
    IF (LTRAJEC) NGCELL = NAIRPARC

    ! SET THE NUMBER OF CELLS PER CPU IN *ATTILA*
    CALL get_dc_index(NGCELL, IDX)
    np = SIZE(IDX,1)
    ALLOCATE(V_NCELL(0:np-1))
    V_NCELL(0:np-1) = IDX(0:np-1,2) - IDX(0:np-1,1) + 1
    NCELL = V_NCELL(p_pe)

    CALL ATTILA_INICOM_2(status)
    IF (status /= 0) CALL error_bi('ATTILA_INICOM_2 reported an error !' &
         , substr)

    LPROC(RANDI_TURB) = LLTBLTURB
    LPROC(RANDI_CONV) = LLCONV
    LPROC(RANDI_CAT)  = LLCAT
    LPROC(RANDI_MOCA) = LLTMOCA

    ! NEW DIMENSIONS
    CALL new_dimension(status, DIMID_NGCELL, 'NCELL', NGCELL)
    CALL channel_halt(substr, status)

    ! NEW REPRESENTATIONS
    CALL new_representation(status, LG_ATTILA, 'LG_ATTILA'    &
         , rank = 1, link = 'x---', dctype = DC_IX            &
         , dimension_ids = (/ DIMID_NGCELL /)                 &
         , ldimlen       = (/ NCELL /)                        &
         , axis = 'N---'                                      &
         )
    CALL channel_halt(substr, status)

    nseg = 1
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    start(:,:) = 1
    cnt(:,:) = 1
    meml(:,:) = 1
    memu(:,:) = 1

    start(:,1) = IDX(p_pe,1)
    cnt(:,1)   = NCELL
    meml(:,1)  = 1
    memu(:,1)  = NCELL

    CALL set_representation_decomp(status, LG_ATTILA &
         , start, cnt, memu, meml, piotype=PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! ----------------------------------------

    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    ! ----------------------------------------

    CALL end_message_bi(modstr, 'INITIALIZATION', substr)

  END SUBROUTINE attila_initialize
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE attila_init_memory

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, nlev, ngpblks
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: LG_ATTILA, SCALAR &
                                         , GP_3D_MID, GP_2D_HORIZONTAL &
                                         , GP_3D_INT
    ! MESSy
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute

    !  COMPLICATED POINTER ARITHMETIC
    !  APOS, NPOS, NPOSM1 ARE USED IN RERUN FILES
    !  IN CASE OF RUNNING ATTILA WITH ECHAM5, APOS,... ARE
    !  POINTERS.
    !  THE PTRS POINT WITH HELP OF SOME AUXILIARY PTRS INTO
    !  THE OUTPUT AND RESTART CHANNELS
    !
    !  -------------------                                    ----------------
    !  | MESSY_ATTILA_E5 |                                    | MESSY_ATTILA |
    !  -------------------                                    ----------------
    !
    !                                         LSPAC(:,:,1,1,1) <= APOS(:,:)
    !
    !  LONZ(:)        => PTRLON(:,1,1,1)   => LSPAC(1,:,:,:,:)
    !  LATZ_MATH(:)   => PTRLAT(:,1,1,1)   => LSPAC(2,:,:,:,:)
    !  etaheight(:)   => PTRETA(:,1,1,1)   => LSPAC(3,:,:,:,:)
    !  pressheight(:) => PTRPRESS(:,1,1,1) => LSPAC(4,:,:,:,:)

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr ='attila_init_memory'
    INTEGER                     :: status
    INTEGER                     :: i
    INTEGER                     :: RND_MP_ATTILA

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    ! OPEN NEW OUTPUT CHANNEL
    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)

    ! op_pj_20141018+
    ! names of attributes must have max. 24 characters
    ! modstr_: attila -> 7 -> 17 available
    CALL new_attribute(status, modstr &
         , modstr//'_I_RANDOM_METHOD', i=I_RANDOM_METHOD)
    CALL channel_halt(substr, status)
    !
    CALL new_attribute(status, modstr &
         , modstr//'_I_RANDOM_PARALLEL', i=I_RANDOM_PARALLEL)
    CALL channel_halt(substr, status)
    !
    CALL new_attribute(status, modstr &
         , modstr//'_I_PBLH_METHOD', i=I_PBLH_METHOD)
    CALL channel_halt(substr, status)
    !
    CALL new_attribute(status, modstr, modstr//'_I_VERT', i=I_VERT)
    CALL channel_halt(substr, status)
    IF (I_VERT == 2) THEN
       CALL new_attribute(status, modstr, modstr//'_PRESS_REF' &
            , r=press_ref)
       CALL channel_halt(substr, status)
    END IF
    !
    IF (LLCONV) THEN
       CALL new_attribute(status, modstr, modstr//'_LLCONV', c='TRUE')
       CALL channel_halt(substr, status)
    ELSE
       CALL new_attribute(status, modstr, modstr//'_LLCONV', c='FALSE')
       CALL channel_halt(substr, status)
    END IF
    CALL channel_halt(substr, status)
    !
    IF (LLTBLTURB) THEN
       CALL new_attribute(status, modstr, modstr//'_LLTBLTURB', c='TRUE')
       CALL channel_halt(substr, status)
    ELSE
       CALL new_attribute(status, modstr, modstr//'_LLTBLTURB', c='FALSE')
       CALL channel_halt(substr, status)
    END IF
    CALL channel_halt(substr, status)
    !
    IF (LTRAJEC) THEN
       CALL new_attribute(status, modstr, modstr//'_LTRAJEC', c='TRUE')
       CALL channel_halt(substr, status)
    ELSE
       CALL new_attribute(status, modstr, modstr//'_LTRAJEC', c='FALSE')
       CALL channel_halt(substr, status)
    END IF
    CALL channel_halt(substr, status)
    !
    ! op_pj_20141018-

    ! ALLOCATE AUXILIAR POINTERS USED DURING INTEGRATON
    ALLOCATE(LSPAC(4,NCELL,1,1,1))
    ALLOCATE(LSPACINT(5,NCELL,1,1,1))
    ALLOCATE(LSPACINTM1(5,NCELL,1,1,1))
    ALLOCATE(LSPACPRESSM1(1,NCELL,1,1,1))

    ! ASSIGNING THE APOS POINTERS
    ptrlon   => LSPAC(1,:,:,:,:)
    ptrlat   => LSPAC(2,:,:,:,:)
    ptreta   => LSPAC(3,:,:,:,:)  ! ETA_HEIGHT
    ptrpress => LSPAC(4,:,:,:,:)  ! PRESS_HEIGHT
    !
    ! TM1
    ptrpressm1 => LSPACPRESSM1(1,:,:,:,:)
    !
    ptrintlon  => LSPACINT(1,:,:,:,:)
    ptrintlat  => LSPACINT(2,:,:,:,:)
    ptrintlev  => LSPACINT(3,:,:,:,:)
    ptrintlon_gp => LSPACINT(4,:,:,:,:)
    ptrintlat_gp => LSPACINT(5,:,:,:,:)
    !
    ptrintlonm1 => LSPACINTM1(1,:,:,:,:)
    ptrintlatm1 => LSPACINTM1(2,:,:,:,:)
    ptrintlevm1 => LSPACINTM1(3,:,:,:,:)
    ptrintlonm1_gp => LSPACINTM1(4,:,:,:,:)
    ptrintlatm1_gp => LSPACINTM1(5,:,:,:,:)

    ! ALLOCATE GLOBAL FIELDS NOT IN CHANNELS
    ! NUMBER OF CELLS PER GRID BOX
    ALLOCATE(NCB_P(SIZE(NCB,1),SIZE(NCB,2),SIZE(NCB,3)))
    ALLOCATE(NCBL_P(SIZE(NCBL,1),SIZE(NCBL,2)))
    NCB_P(:,:,:)  = 0.0_dp
    NCBL_P(:,:)   = 0.0_dp
    !
    ALLOCATE(zaps_scb(nproma, ngpblks))
    ALLOCATE(zpi_scb(nproma, nlev+1, ngpblks))

    CALL new_channel_object(status, modstr, 'PLAT-MATH' &
         , p1=latz_math, reprid=LG_ATTILA, mem=ptrlat, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PLAT-MATH' &
         , 'long_name', c='latitude position' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PLAT-MATH' &
         , 'units', c='degreees; 0 (NP) - 180 (SP)' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'PLAT' &
         , p1=latz, reprid=LG_ATTILA, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PLAT' &
         , 'long_name', c='latitude position' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PLAT' &
         , 'units', c='degrees_north' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'PLON' &
         , p1=lonz, reprid=LG_ATTILA, mem=ptrlon, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PLON' &
         , 'long_name', c='longitude position' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PLON' &
         , 'units', c='degrees_east' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'PETA' &
         , p1=etaheight, reprid=LG_ATTILA, mem=ptreta, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PETA' &
         , 'long_name', c='level position (eta coord.)' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PETA' &
         , 'units', c='1' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'PPRESS' &
         , p1=pressheight, reprid=LG_ATTILA, mem=ptrpress, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PPRESS' &
         , 'long_name', c='level position (pressure)' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PPRESS' &
         , 'units', c='Pa' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'PPRESS_M1' &
         , p1=press_m1, reprid=LG_ATTILA, mem=ptrpressm1, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PPRESS_M1' &
         , 'long_name', c='level position (pressure) at t-1' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PPRESS_M1' &
         , 'units', c='Pa' )
    CALL channel_halt(substr, status)

    ! INTEGER POSITION VALUES
    CALL new_channel_object(status, modstr, 'IPLON' &
         , reprid=LG_ATTILA, mem=ptrintlon_gp, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'IPLON' &
         , 'long_name', c='longitude position (index)' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'IPLON' &
         , 'units', c=' ' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'IPLAT' &
         , reprid=LG_ATTILA, mem=ptrintlat_gp, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'IPLAT' &
         , 'long_name', c='latitude position (index)' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'IPLAT' &
         , 'units', c=' ' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'IPLEV' &
         , reprid=LG_ATTILA, mem=ptrintlev, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'IPLEV' &
         , 'long_name', c='level position (index)' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'IPLEV' &
         , 'units', c=' ' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'IPLON_A' &
         , reprid=LG_ATTILA, mem=ptrintlon, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'IPLAT_A' &
         , reprid=LG_ATTILA, mem=ptrintlat, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)

    ! INTEGER POSITION VALUES AT t=tm1
    CALL new_channel_object(status, modstr, 'IPLON_M1' &
         , reprid=LG_ATTILA, mem=ptrintlonm1_gp, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'IPLON_M1' &
         , 'long_name', c='longitude position (index) at t-1' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'IPLON_M1' &
         , 'units', c=' ' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'IPLAT_M1' &
         , reprid=LG_ATTILA, mem=ptrintlatm1_gp, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'IPLAT_M1' &
         , 'long_name', c='latitude position (index) at t-1' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'IPLAT_M1' &
         , 'units', c=' ' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'IPLEV_M1' &
         , reprid=LG_ATTILA, mem=ptrintlevm1, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'IPLEV_M1' &
         , 'long_name', c='level position (index) at t-1' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'IPLEV_M1' &
         , 'units', c=' ' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'IPLON_A_M1' &
         , reprid=LG_ATTILA, mem=ptrintlonm1, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'IPLAT_A_M1' &
         , reprid=LG_ATTILA, mem=ptrintlatm1, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)

    ! ALLOCATION AND POINTER ASSIGNMENTS IN *ATTILA*
    CALL ATTILA_ALLOC(status, &
         LSPAC=LSPAC,LSPACINT=LSPACINT,LSPACINTM1=LSPACINTM1, &
         LSPACPRESSM1=LSPACPRESSM1)
    IF (status /= 0) &
         CALL error_bi('MEMORY ALLOCATION/POINTER ASSOCIATION FAILED!', substr)

    !
    ! ADDING SOME MORE IMPORTANT VARIABLES TO THE CHANNEL
    !
    ! NUMBER OF CELLS PER GRID BOX
    CALL new_channel_object(status, modstr, 'NCB' &
         , p3=SGNCB, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NCB' &
         , 'long_name', c='number of cells in each grid box' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NCB' &
         , 'units', c=' ' )
    CALL channel_halt(substr, status)
    ! NUMBER OF CELLS PER GRID BOX AT T-1
    CALL new_channel_object(status, modstr, 'NCBM1' &
         , p3=SGNCBM1, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NCBM1' &
         , 'long_name', c='number of cells in each grid box' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NCBM1' &
         , 'units', c=' ' )
    CALL channel_halt(substr, status)
    ! NUMBER OF CELLS PER GRID BOX IN BOUNDARY LAYER
    CALL new_channel_object(status, modstr, 'NCBL' &
         , p2=SGNCBL, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NCBL' &
         , 'long_name', c='number of cells in BL' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NCBL' &
         , 'units', c=' ' )
    CALL channel_halt(substr, status)

    ! MASS OF ATTILA CELLS
    CALL new_channel_object(status, modstr, 'AMCELL' &
         , p0=ptramcell, reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'AMCELL' &
         , 'long_name', c='mass of CELL' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'AMCELL' &
         , 'units', c='kg' )
    CALL channel_halt(substr, status)

    ! op_sb_20140211+
    IF (I_VERT == 2) THEN
       CALL new_channel_object(status, modstr, 'tempdot', p3=tempdot &
            , reprid=GP_3D_MID )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'tempdot', 'long_name', c='tendency diabatic processes' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'tempdot', 'units', c='K/s' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'dthetadt', p3=dthetadt &
            , reprid=GP_3D_MID )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'dthetadt', 'long_name', c='tendency (theta) diabatic processes' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'dthetadt', 'units', c='K/s' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'xidot', p3=xidot &
            , reprid=GP_3D_INT )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'xidot', 'long_name', c='vertical velocity in isentropes' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'xidot', 'units', c='K/s')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'xidot_cor', p3=xidot_cor &
            , reprid=GP_3D_INT )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'xidot_cor', 'long_name', c='vertical velocity in isentropes cfl-criterion corrected' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'xidot_cor', 'units', c='K/s')
       CALL channel_halt(substr, status)

       ! op_sb_20160707+
       CALL new_channel_object(status, modstr, 'sidot', p3=sidot &
            , reprid=GP_3D_INT )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'sidot', 'long_name', c='vertical velocity in sigma' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'sidot', 'units', c='')
       CALL channel_halt(substr, status)
       ! op_sb_20160707-

       CALL new_channel_object(status, modstr, 'diab', p3=diab &
            , reprid=GP_3D_MID )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'diab', 'long_name', c='' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'diab', 'units', c='' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'kine', p3=kine &
            , reprid=GP_3D_MID )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'kine', 'long_name', c='' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'kine', 'units', c='K/s' )
       CALL channel_halt(substr, status)

     END IF
    ! op_sb_20140211-

    ! FLAG FOR CONVECTION ALONG TRAJECTORY
    IF (LLCONV) THEN
       CALL new_channel_object(status, modstr, 'CONVDELK' &
            , p1=CONVDELK, reprid=LG_ATTILA)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'CONVDELK' &
            , 'long_name', c='level index change by convection' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'CONVDELK' &
            , 'units', c=' ' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'NCONV_ITER' &
            , p0=NCONV_ITER, reprid=SCALAR)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'NCONV_ITER' &
            , 'long_name', c='number of convection iterations' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'NCONV_ITER' &
            , 'units', c=' ' )
       CALL channel_halt(substr, status)
       !
       CALL new_channel_object(status, modstr, 'NCB_CUD' &
            , p3=SGNCB_CUD, reprid=GP_3D_MID, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'NCB_CUD' &
            , 'long_name', c='total number of parcels after up/downdraft' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'NCB_CUD' &
            , 'units', c=' ' )
       CALL channel_halt(substr, status)
       !
       CALL new_channel_object(status, modstr, 'NCB_CSUB' &
            , p3=SGNCB_CSUB, reprid=GP_3D_MID, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'NCB_CSUB' &
            , 'long_name', c='total number of parcels after subsidence' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'NCB_CSUB' &
            , 'units', c=' ' )
       CALL channel_halt(substr, status)
       !
       CALL new_channel_object(status, modstr, 'pos_beg_up' &
            , p1=pos_beg_up, reprid=LG_ATTILA)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'pos_beg_up' &
            , 'long_name', c='level of parcel before updraft' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'pos_beg_up' &
            , 'units', c=' ' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'pos_end_up' &
            , p1=pos_end_up, reprid=LG_ATTILA)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'pos_end_up' &
            , 'long_name', c='level of parcel after updraft' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'pos_end_up' &
            , 'units', c=' ' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'pos_end_do' &
            , p1=pos_end_do, reprid=LG_ATTILA)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'pos_end_do' &
            , 'long_name', c='level of parcel after downdraft' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'pos_end_do' &
            , 'units', c=' ' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'pos_end_sub' &
            , p1=pos_end_sub, reprid=LG_ATTILA)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'pos_end_sub' &
            , 'long_name', c='level of parcel after subsidence' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'pos_end_sub' &
            , 'units', c=' ' )
       CALL channel_halt(substr, status)

       ! MASS OF ATTILA CELLS
       CALL new_channel_object(status, modstr, 'rlg_conv' &
            , p0=rlg_conv, reprid=SCALAR, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'rlg_conv' &
            , 'long_name', c='convection flag' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'rlg_conv' &
            , 'units', c='' )
       CALL channel_halt(substr, status)
       rlg_conv = 1._dp

    END IF

    IF (LVDIAG) THEN
       ! ADDITIONAL VELOCITY DIAGNOSTICS
       CALL new_channel_object(status, modstr, 'uvel' &
            , p1=uvel, reprid=LG_ATTILA)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'uvel' &
            , 'long_name', c='u-velocity' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'uvel' &
            , 'units', c='m/s' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'vvel' &
            , p1=vvel, reprid=LG_ATTILA)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'vvel' &
            , 'long_name', c='v-velocity' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'vvel' &
            , 'units', c='m/s' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'wvel' &
            , p1=wvel, reprid=LG_ATTILA)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'wvel' &
            , 'long_name', c='w-velocity' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'wvel' &
            , 'units', c='1/s' )
       CALL channel_halt(substr, status)
    END IF

    ! SPECIAL IN TRAJECTORY MODE
    IF (LTRAJEC) THEN
       CALL new_channel_object(status, modstr, 'RTRAJEC_INTEGRATE' &
            , p1=RTRAJEC_INTEGRATE, reprid=LG_ATTILA, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'RTRAJEC_INTEGRATE' &
            , 'long_name', c='time integration flag' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'RTRAJEC_INTEGRATE' &
            , 'units', c=' ' )
       CALL channel_halt(substr, status)
       ! INTITIALIZE
       RTRAJEC_INTEGRATE(:) = -1.0_dp
    END IF

    IF (L_RANDOM_TEST) THEN
       CALL new_channel_object(status, modstr, 'test_rnd' &
            , p1=test_rnd, reprid=LG_ATTILA)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'test_rnd' &
            , 'long_name', c='test of random number sequence' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'test_rnd' &
            , 'units', c=' ' )
       CALL channel_halt(substr, status)
    END IF

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    CALL start_message_bi(modstr, 'RANDOM NUMBER INITIALISATION', substr)

    ! ALLOCATE SPACE FOR RANDOM NUMBERS AND
    ! INITIALISE PSEUDO RANDOM NUMBER STREAMS.
    !
    ! NOTE: DO THIS ALWAYS FOR TURB, MC, CONV AND CAT, EVEN IF
    !       LLTBLTURB, LLTMOCA, LLCONV AND/OR LLCAT ARE .FALSE.;
    !       THIS ALLOWS A MORE SIMPLER CODE FOR INTEGRATION AND
    !       SWITCHING ON/OFF PROCESSES AFTER RESTART WITHOUT
    !       MESSING UP THE DIFFERENT RANDOM NUMBER STREAMS.
    !
    ! ALLOCATE WORKSPACE FOR RANDOM NUMBERS
    DO i=1, NRPROC
!qqq TEMPORARY SPECIAL CASE FOR CONV
       IF (i == RANDI_CONV) THEN
          ALLOCATE(HARVEST(i)%ptr(NGCELL,RAND_MULT(i)))
          RND_MP_ATTILA = RND_MP_PIN
       ELSE
          ALLOCATE(HARVEST(i)%ptr(NCELL,RAND_MULT(i)))
          RND_MP_ATTILA = I_RANDOM_PARALLEL
       ENDIF

       ! use default jump of p_pe*2^256 for RND_ATTILA = RND_MP_PSJ
       CALL rnd_init_bi(rndid(i), RND_ATTILA, RND_MP_ATTILA, START_SEED(i))
    END DO

    CALL end_message_bi(modstr, 'RANDOM NUMBER INITIALISATION', substr)

  END SUBROUTINE attila_init_memory
!---------------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE attila_init_coupling(flag)

    ! ATTILA_E5 MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! initialize 'coupling' to online tracers/channels
    !
    ! Author: Michael Traub, MPICH, Feb 2004

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: LG_FMM_POS
    ! MESSy
    USE messy_main_channel,          ONLY: get_channel_object &
                                         , new_channel_object &
                                         , new_attribute
#ifdef MESSYTENDENCY
    USE messy_main_tendency_bi,      ONLY: mtend_request, mtend_id_t
#endif

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr = 'attila_init_coupling'
    INTEGER                           :: status


    SELECT CASE(flag)

    CASE(1)
    CALL start_message_bi(modstr,'COUPLING INITIALIZATION',substr)

    IF (L_INI_PARALLEL) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) 'Parallel initialisation ...'
       ENDIF
    ELSE
       IF (p_parallel_io) THEN
          WRITE(*,*) 'Single PE initialisation ...'
       ENDIF
    END IF

    ! TURBULENCE IN PLANETARY BOUNDARY LAYER
    IF (LLTBLTURB) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) 'Checking for planetary boundary layer index ... '
          WRITE(*,*) '    channel: ',TRIM(C_PBLH_INDEX(1))
          WRITE(*,*) '    object : ',TRIM(C_PBLH_INDEX(2))
       END IF
       CALL get_channel_object(status &
            , TRIM(C_PBLH_INDEX(1)), TRIM(C_PBLH_INDEX(2)) &
            , p2=pblh_i )
       CALL channel_halt(substr, status)
    END IF

    ! CLEAR AIR TURBULENCE INDEX
    IF (LLCAT) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) 'Checking for clear air turbulence index ... '
          WRITE(*,*) '    channel: ',TRIM(C_TURB_CATI(1))
          WRITE(*,*) '    object : ',TRIM(C_TURB_CATI(2))
       END IF
       CALL get_channel_object(status &
            , TRIM(C_TURB_CATI(1)), TRIM(C_TURB_CATI(2)) &
            , p3=cati )
       CALL channel_halt(substr, status)
    END IF

    ! COUPLING TO GP CONVECTION PARAMETERS
    IF (LLCONV) THEN

       ! 1) UPDRAFT MASS FLUX
       IF (p_parallel_io) THEN
          WRITE(*,*) 'Checking for convective updraft mass flux ... '
          WRITE(*,*) '    channel: ',TRIM(C_CONV_uflx(1))
          WRITE(*,*) '    object : ',TRIM(C_CONV_uflx(2))
       END IF
       CALL get_channel_object(status &
            , TRIM(C_CONV_uflx(1)), TRIM(C_CONV_uflx(2)) &
            , p3=uflx_conv )
       CALL channel_halt(substr, status)

       ! 2) UPDRAFT MASS FLUX
       IF (p_parallel_io) THEN
          WRITE(*,*) 'Checking for convective downdraft mass flux ... '
          WRITE(*,*) '    channel: ',TRIM(C_CONV_dflx(1))
          WRITE(*,*) '    object : ',TRIM(C_CONV_dflx(2))
       END IF
       CALL get_channel_object(status &
            , TRIM(C_CONV_dflx(1)), TRIM(C_CONV_dflx(2)) &
            , p3=dflx_conv )
       CALL channel_halt(substr, status)

       ! 3) TYPE OF CONVECTION
       IF (p_parallel_io) THEN
          WRITE(*,*) 'Checking for type of convection ... '
          WRITE(*,*) '    channel: ',TRIM(C_CONV_type(1))
          WRITE(*,*) '    object : ',TRIM(C_CONV_type(2))
       END IF
       CALL get_channel_object(status &
            , TRIM(C_CONV_type(1)), TRIM(C_CONV_type(2)) &
            , p2=type_conv )
       CALL channel_halt(substr, status)
    END IF

    ! op_sb_20130823+
    IF (I_VERT /= 1) THEN
       CALL get_channel_object(status, 'ECHAM5', 'sigmadot', p3=sigmadot)
       CALL channel_halt(substr, status)
    END IF

    IF (I_VERT == 2) THEN
! op_sb_20160923+
! Due to large short term variability, the on-line diagnosed tropopause
! cannot be used here, because it causes air parcel overshoots.
!!$    CALL get_channel_object(status, 'tropop', 'tp', p2=tp)
       CALL get_channel_object(status, 'tropop', 'tp_clim', p2=tp)
! op_sb_20160923-
       CALL channel_halt(substr, status)
!!$#ifdef MESSYTENDENCY
!!$       CALL mtend_request('dyn', mtend_id_t, tte_dyn)
!!$#endif
    ENDIF
    ! op_sb_20130823-

    ! these objects are only required, if LGFMM is switched on
    CALL get_channel_object(status, 'lgfmm', 'slgmass')
    L_LG_FMM = (status == 0)
    IF (L_LG_FMM) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) 'Creating additional channel object for lgfmm ...'
       END IF
       CALL new_channel_object(status, modstr, 'APOS_FMM' &
            , p3=APOS_FMM, reprid=LG_FMM_POS)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'APOS_FMM' &
            , 'long_name', c='position information (RK4); NCELL x [xyz] x [T]' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'APOS_FMM' &
            , 'units', c='degreas_east ; degrees_north ; eta ' )
     ENDIF

    CALL end_message_bi(modstr,'COUPLING INITIALIZATION',substr)
 CASE(2)

#ifdef MESSYTENDENCY
    IF (I_VERT == 2) THEN
       CALL mtend_request('dyn', mtend_id_t, tte_dyn)
    END IF
#endif

    CASE DEFAULT
       CALL error_bi('UNKNOWN FLAG', substr)

    END SELECT
  END SUBROUTINE attila_init_coupling
! ----------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE attila_init_tracer

    IMPLICIT NONE

    ! INITIALIZE POSITIONS HERE, BECAUSE IT IS NEEDED FOR TRACER
    ! INITIALIZATION
    IF (.NOT.LTRAJEC) CALL attila_initialize_positions

  END SUBROUTINE attila_init_tracer
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE attila_local_start

    ! ECHAM5/MESSy
    USE messy_main_data_bi,         ONLY: tte_3d
    USE messy_main_grid_def_mem_bi, ONLY: jrow

    IMPLICIT NONE

#ifndef MESSYTENDENCY
    IF (I_VERT /= 2) RETURN
    ! save temperature tendency after advection
    tempdot(:,:,jrow) = tte_3d(:,:,jrow)
#endif

  END SUBROUTINE attila_local_start
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE attila_local_end

    ! ECHAM5/MESSy
    USE messy_main_data_bi,         ONLY: tte_3d
    USE messy_main_grid_def_mem_bi, ONLY: jrow

    IMPLICIT NONE

    IF (I_VERT /= 2) RETURN
    ! calculate dtdt for diabatic processes only
#ifndef MESSYTENDENCY
    tempdot(:,:,jrow) = tte_3d(:,:,jrow) - tempdot(:,:,jrow)
#else
    tempdot(:,:,jrow) = tte_3d(:,:,jrow) - tte_dyn(:,:,jrow)
#endif

  END SUBROUTINE attila_local_end
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
! op_sb_20140516+
! ... everything moved to global_end
!!$  SUBROUTINE attila_global_start
!!$  END SUBROUTINE attila_global_start
! op_sb_20140516+
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE attila_global_end

    ! ECHAM5/MESSy
    USE messy_main_transform_bi,    ONLY: trp_gpdc_gpgl, get_dc_index
    USE messy_main_data_bi,         ONLY: etadot_3d       &   ! scan1
                                        , tsurf_2d        &
                                        , rinum_3d        &
                                        , u_scb, v_scb    &   ! scan_buffer
                                        , t_scb, alps_scb &   ! scan_buffer
                                        , rh_scb          &   ! scan_buffer
                                        , tpot_3d
    USE messy_main_grid_def_bi,     ONLY:sqcst_2d
    USE messy_main_grid_def_mem_bi, ONLY: vct, nlevp1, nvclev, nlev &
                                        , nproma, ngpblks
    USE messy_main_timer,           ONLY: lstart, lresume
    ! MESSSy
    USE messy_main_constants_mem,   ONLY: cp_air, R_gas, M_air

    IMPLICIT NONE

    REAL(DP), PARAMETER :: rd = 1000._DP * R_gas / M_air ! J/mol/kg

    INTEGER :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'attila_global_end'
    REAL(dp), DIMENSION(:,:,:), POINTER :: zu_scb  => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: zv_scb  => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: zt_scb  => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: zrh_scb => NULL()
    REAL(dp), DIMENSION(:,:)  , POINTER :: tsurf => NULL() ! surface temperature
    REAL(dp), DIMENSION(:,:,:), POINTER :: rinr  => NULL() ! Richardson number
    INTEGER                             :: jk
    REAL(dp), DIMENSION(NLON,NGL)       :: zrcst

    L_INI_FINI = .NOT. lstart ! op_sb_20141015

    zu_scb => u_scb
    zv_scb => v_scb
    zt_scb => t_scb
    zrh_scb => rh_scb

    ! SURFACE PRESSURE
!qqq
    ! WHAT ABOUT surface pressure tendency ???
    zaps_scb(:,:)  = exp(alps_scb(:,:))          ! [Pa]
    ! PRESSURE AT LAYER INTERFACES
    DO jk=1, nlevp1
       zpi_scb(:,jk,:) = vct(jk) + vct(nvclev+jk) * zaps_scb(:,:) ! [Pa]
    END DO

    ! ALWAYS NEEDED
    CALL trp_gpdc_gpgl(1, zpi_scb,   pph)
    CALL trp_gpdc_gpgl(1, zu_scb,    pwu)
    CALL trp_gpdc_gpgl(1, zv_scb,    pwv)
    CALL trp_gpdc_gpgl(1, zaps_scb,  pps)
    CALL trp_gpdc_gpgl(1, sqcst_2d,  sqcst)

    ! op_sb_20140121+
    SELECT CASE (I_VERT)
    CASE(1)  ! default
       CALL trp_gpdc_gpgl(1, etadot_3d, pww)
    CASE(2)  ! thetadot hybrid
       CALL trp_gpdc_gpgl(1, zt_scb,    ptemp)
       CALL trp_gpdc_gpgl(1, tpot_3d,   ptpot) ! op_sb_20140211

       CALL attila_calc_sidot                  ! op_sb_20160707
       CALL attila_calc_xidot                  ! op_sb_20140211

       CALL trp_gpdc_gpgl(1, xidot,     pww)
       CALL trp_gpdc_gpgl(1, sigma_ref,  sigma_ref_g) ! op_sb_20151113
       CALL trp_gpdc_gpgl(1, i_ref,  i_ref_g)
       ! sb_op_20160503+
       IF (lstart) THEN
          tsurf => zt_scb(:,nlev,:)
         ELSE
          tsurf => tsurf_2d
       ENDIF
       CALL trp_gpdc_gpgl(1, tsurf,     surftemp)
       ! sb_op_20160503-

    CASE(3)  ! sigmadot
       CALL trp_gpdc_gpgl(1, sigmadot, pww)
    END SELECT
    ! op_sb_20140121-

    ! DIVIDE (U,V)*cos(lat) BY cos(lat) TO GET HORIZONTAL WIND
    zrcst(:,:) = 1._dp/ sqcst(:,:)
    DO jk = 1, nlev
       pwu(:,jk,:) = pwu(:,jk,:) * zrcst(:,:)
       pwv(:,jk,:) = pwv(:,jk,:) * zrcst(:,:)
    END DO

    ! FIRST LAYER IN FREE ATMOSPHERE
    SELECT CASE (I_PBLH_METHOD)
    CASE(0) ! ------------- INTERNAL CALCULATION ------------
       !
       ! RICHARDSON NUMBER AND SURFACE TEMPERATURE
       ! quick and dirty workaround for 1st time step
       IF (lstart) THEN
          tsurf => zt_scb(:,nlev,:)
          ALLOCATE(rinr(nproma, nlev, ngpblks))
          rinr(:,:,:) = 1.0_dp ! ???
       ELSE
          tsurf => tsurf_2d
          rinr  => rinum_3d
       END IF
       !
       CALL trp_gpdc_gpgl(1, tsurf,     surftemp)
       CALL trp_gpdc_gpgl(1, zt_scb,    ptemp)    ! possibly called twice !
       !                                          ! (see above)
       CALL trp_gpdc_gpgl(1, zrh_scb,   geop_3d)
       CALL trp_gpdc_gpgl(1, tpot_3d,   pteta1_3d) ! possibly called twice !
       !                                           ! (see above)
       CALL trp_gpdc_gpgl(1, rinr,      priad_3d)
       ! CLEAN MEMORY
       IF (lstart) THEN
          IF (ASSOCIATED(rinr)) DEALLOCATE(rinr)
          NULLIFY(rinr)
       END IF
       !
       ! NOTE: THIS ROUTINE HAS SOME PROBLEMS AND IS NOT RECOMMENDET
       CALL ATTILA_FIRSTL(surftemp, geop_3d, pteta1_3d, priad_3d, tpot_3d)
       !
    CASE(1) ! ------------- EXTERNAL CALCULATION ------------
       !
       IF (LLTBLTURB) THEN
          IF (.NOT.lstart) THEN
             !
             !NOTE: pblh_i in TROPOP is calculated in tropo_vdiff!
             !      Therefore, calling ATTILA_FIRSTL_EXT at lstart would cause
             !      resetting of KHPBL (messy_attila.f90) to ZERO instead
             !      of nlev-2 ...
             !
             CALL trp_gpdc_gpgl(1, pblh_i, pblh_idx_gl)
             CALL ATTILA_FIRSTL_EXT(status, pblh_idx_gl)
             IF (status /= 0) &
                  CALL error_bi('ATTILA_FIRSTL_EXT reported an error',substr)
          END IF
       END IF
       !
    CASE DEFAULT ! -------- ERROR ---------------------------
       !
       CALL error_bi(' *** ERRROR *** UNKNOWN METHOD OF PBLH-CALC.', substr)
       !
    END SELECT

    ! SETUP FIELDS
    CALL attila_fill(status)
    IF (status /= 0) &
         CALL error_bi('ATTILA_FILL reported an error', substr)

    IF (LTRAJEC) CALL attila_initialize_positions_t(1) ! op_pj_20141020

    IF (lstart .or. lresume) THEN ! op_sb_20141031 (.OR. lresume)
       ! UPDATE CELL DISTRIBUTION
       CALL attila_update_celldist
       L_INI_FINI = .TRUE.            ! op_sb_20140930
    END IF

    SGNCBM1(:,:,:) = SGNCB(:,:,:)  ! save value for ltm1 trafo

    CALL attila_integrate

!------original global_start ends here   ! op_sb_20140516-

    IF (LLCONV) THEN
       CALL attila_convect
       CALL attila_update_celldist
    END IF

    IF (LTRAJEC) CALL attila_initialize_positions_t(2)

    IF (I_VERT == 2) CALL trp_gpdc_gpgl(-1, xidot_cor,     pww)

 ! op_sb_20151113+
   IF (ASSOCIATED(sigma_ref)) THEN
       DEALLOCATE(sigma_ref) ; NULLIFY(sigma_ref)
   ENDIF
   IF (ASSOCIATED(i_ref)) THEN
       DEALLOCATE(i_ref) ; NULLIFY(i_ref)
   ENDIF
 ! op_sb_20151113-

  END SUBROUTINE attila_global_end
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE attila_free_memory

    IMPLICIT NONE

    ! LOCAL
    INTEGER :: i

    CALL attila_global_exit  ! 'INVERSE' OF SUBROUTINE ATTILA_ALLOC

    IF (ASSOCIATED(LSPAC)) THEN
       DEALLOCATE(LSPAC) ; NULLIFY(LSPAC)
    END IF
    IF (ASSOCIATED(LSPACINT)) THEN
       DEALLOCATE(LSPACINT) ; NULLIFY(LSPACINT)
    END IF
    IF (ASSOCIATED(LSPACINTM1)) THEN
       DEALLOCATE(LSPACINTM1) ; NULLIFY(LSPACINTM1)
    END IF
    IF (ASSOCIATED(LSPACPRESSM1)) THEN
       DEALLOCATE(LSPACPRESSM1) ; NULLIFY(LSPACPRESSM1)
    END IF

    IF (ASSOCIATED(NCB_P)) THEN
       DEALLOCATE(NCB_P) ; NULLIFY(NCB_P)
    END IF
    IF (ASSOCIATED(NCBL_P)) THEN
       DEALLOCATE(NCBL_P) ; NULLIFY(NCBL_P)
    END IF
    IF (ASSOCIATED(IDX)) THEN
       DEALLOCATE(IDX) ; NULLIFY(IDX)
    END IF
    IF (ASSOCIATED(V_NCELL)) THEN
       DEALLOCATE(V_NCELL) ; NULLIFY(V_NCELL)
    END IF

    ! RANDOM NUMBERS AND INDICES
    DO i=1, NRPROC
       IF (ASSOCIATED(HARVEST(i)%ptr)) DEALLOCATE(HARVEST(i)%ptr)
       NULLIFY(HARVEST(i)%ptr)
       CALL rnd_finish_bi(rndid(i))
    END DO

    ! LOCAL FIELDS
    IF (ASSOCIATED(zaps_scb)) THEN
       DEALLOCATE(zaps_scb) ; NULLIFY(zaps_scb)
    END IF
    IF (ASSOCIATED(zpi_scb)) THEN
       DEALLOCATE(zpi_scb) ; NULLIFY(zpi_scb)
    END IF

    ! DELETE GLOBALIZED FIELDS ON ANY CPU
    IF (ASSOCIATED(pwu)) THEN
       DEALLOCATE(pwu) ; NULLIFY(pwu)
    ENDIF
    IF (ASSOCIATED(pwv)) THEN
       DEALLOCATE(pwv) ; NULLIFY(pwv)
    ENDIF
    IF (ASSOCIATED(pph)) THEN
       DEALLOCATE(pph) ; NULLIFY(pph)
    ENDIF
    IF (ASSOCIATED(pww)) THEN
       DEALLOCATE(pww) ; NULLIFY(pww)
    ENDIF
    IF (ASSOCIATED(sqcst)) THEN
       DEALLOCATE(sqcst) ; NULLIFY(sqcst)
    END IF
    IF (ASSOCIATED(ptpot)) THEN
       DEALLOCATE(ptpot) ; NULLIFY(ptpot)
    ENDIF
    ! op_sb_20141113+
    IF (ASSOCIATED(sigma_ref_g)) THEN
       DEALLOCATE(sigma_ref_g) ; NULLIFY(sigma_ref_g)
    ENDIF
    ! op_sb_20141113-
    ! op_sb_20150710+
    IF (ASSOCIATED(i_ref_g)) THEN
       DEALLOCATE(i_ref_g) ; NULLIFY(i_ref_g)
    ENDIF
    ! op_sb_20160710-

    ! -> INTERNAL CALCULATION OF PBLH
    IF (ASSOCIATED(surftemp)) THEN
       DEALLOCATE(surftemp) ; NULLIFY(surftemp)
    END IF
    IF (ASSOCIATED(ptemp)) THEN
       DEALLOCATE(ptemp) ; NULLIFY(ptemp)
    ENDIF

    IF (ASSOCIATED(geop_3d)) THEN
       DEALLOCATE(geop_3d) ; NULLIFY(geop_3d)
    ENDIF
    IF (ASSOCIATED(pteta1_3d)) THEN
       DEALLOCATE(pteta1_3d) ; NULLIFY(pteta1_3d)
    ENDIF
    IF (ASSOCIATED(priad_3d)) THEN
       DEALLOCATE(priad_3d) ; NULLIFY(priad_3d)
    ENDIF
    ! -> EXTERNAL CALCULATION OF PBLH
    IF (ASSOCIATED(pblh_idx_gl)) THEN
       DEALLOCATE(pblh_idx_gl) ; NULLIFY(pblh_idx_gl)
    ENDIF
    ! -> CONVECTION
    IF (ASSOCIATED(umassfl)) THEN
       DEALLOCATE(umassfl) ; NULLIFY(umassfl)
    ENDIF
    IF (ASSOCIATED(dmassfl)) THEN
       DEALLOCATE(dmassfl) ; NULLIFY(dmassfl)
    ENDIF
    IF (ASSOCIATED(typeconv)) THEN
       DEALLOCATE(typeconv) ; NULLIFY(typeconv)
    ENDIF

  END SUBROUTINE attila_free_memory
!---------------------------------------------------------------------------

! ===========================================================================
! PRIVATE ATTILA INTERFACE ROUTINES
! ===========================================================================

! ----------------------------------------------------------------------
  SUBROUTINE attila_read_nml_cpl(status, iou)

    ! ATTILA MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5/ATTILA
    !
    ! Author: Patrick Joeckel, MPICH, Oct 2003

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'attila_read_nml_cpl'

    NAMELIST /CPL/ &
         L_INI_PARALLEL, I_RANDOM_METHOD, I_RANDOM_PARALLEL, L_RANDOM_TEST &
         , C_PBLH_INDEX  &
         , C_CONV_UFLX,  C_CONV_DFLX, C_CONV_ENTR, C_CONV_DNTR,  C_CONV_TPCV &
         , C_CONV_TYPE, C_TURB_CATI

    NAMELIST /TRAJ/ AP

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status
    INTEGER              :: jn

    ! ERROR STATUS
    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)

    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST

    CALL read_nml_close(substr, iou, modstr)

    ! TRAJECTORY MODE
    IF (LTRAJEC) THEN

       DO jn=1,MAXNAIRPARC
          AP(jn) = POSITION( 0.0_dp, 0.0_dp, -1.0_dp, 0, 0, 0, 0, 0 )
       END DO

       CALL read_nml_open(lex, substr, iou, 'TRAJ', modstr)
       IF (.not.lex) RETURN    ! <modstr>.nml does not exist

       READ(iou, NML=TRAJ, IOSTAT=fstat)

       CALL read_nml_check(fstat, substr, iou, 'TRAJ', modstr)
       IF (fstat /= 0) RETURN  ! error while reading namelist

       DO jn = 1, MAXNAIRPARC
!          WRITE(*,*) jn, AP(jn)
          IF (AP(jn)%press < 0.0_dp) EXIT
       END DO
       NAIRPARC = jn -1

       WRITE(*,*) 'NUMBER OF INDIVIDUAL TRAJECTORIES: ',NAIRPARC

       ! SET START TIME TO DEFAULT, IF REQUESTED
       IF (LTRAJEC_SAME_DATE) THEN
          DO jn = 1, NAIRPARC
             AP(jn)%year  = LTRAJEC_DATE(1)
             AP(jn)%month = LTRAJEC_DATE(2)
             AP(jn)%day   = LTRAJEC_DATE(3)
             AP(jn)%hour  = LTRAJEC_DATE(4)
             AP(jn)%minute = LTRAJEC_DATE(5)
          END DO
       END IF

       CALL read_nml_close(substr, iou, modstr)
    END IF

    status = 0  ! no ERROR

  END SUBROUTINE attila_read_nml_cpl
! ----------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE attila_initialize_positions

    ! called in attila_init_tracer

    ! ECHAM5/MESSy
    USE messy_main_transform_bi, ONLY: scatter_glix
    USE messy_main_mpi_bi,       ONLY: p_io, p_pe
    USE messy_main_timer,        ONLY: lresume

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'attila_initialize_positions'
    INTEGER                     :: status
    INTEGER                     :: i, j
    ! START SEED VALUES (INITIALIZATION ONLY)
    INTEGER, PARAMETER          :: ISEED1 = 123321
    INTEGER, PARAMETER          :: ISEED2 = 8472
    ! RANDOM NUMBER ARRAYS (REQUIRED ONCE PER SIMULATION FOR INITIALISATION)
    REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: HARVESTINIPOS1 ! ATTILA_INIPOS
    REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: HARVESTINIPOS2 ! FOR ADJUSTMENT
    INTEGER                                :: rndid1, rndid2 ! rnd stream ids
    !
    INTEGER :: NCELL_SAVE
    REAL(dp), DIMENSION(:,:), POINTER, SAVE :: TMP_APOS   => NULL() ! 4xNCELL
    REAL(dp), DIMENSION(:,:), POINTER, SAVE :: TMP_NPOS   => NULL() ! 5xNCELL
    REAL(dp), DIMENSION(:,:), POINTER, SAVE :: TMP_NPOSM1 => NULL() ! 5xNCELL
    REAL(dp), DIMENSION(:),   POINTER :: tmp_l => NULL()
    REAL(dp), DIMENSION(:),   POINTER :: tmp_g => NULL()

    ! NO INITIALISATION REQUIRED AFTER RESTART
    init: IF (.NOT.lresume) THEN

       ! ### RADNOM POSITIONING AT VERY FIRST TIME STEP #################

       l_ini_par: IF (L_INI_PARALLEL) THEN

          !
          ! NOTE: This block contains the ('old') parallel
          !       initialisation, which has the drawback that it is dependent
          !       on the chosen number of CPUs and results in inhomogeneous
          !       initial cell distributions for #CPUs > 1.
          !
          CALL rnd_init_bi(rndid1, RND_ATTILA, RND_MP_PSJ, ISEED1)
          CALL rnd_init_bi(rndid2, RND_ATTILA, RND_MP_PSJ, ISEED2)

          ALLOCATE(HARVESTINIPOS1(NCELL,NR1))
          ALLOCATE(HARVESTINIPOS2(NCELL,NR2))

          DO j=1, NR1
             CALL rnd_number_bi(rndid1, HARVESTINIPOS1(:,j) &
                  , ng=NGCELL, nl=V_NCELL(0:np-1) )
          END DO

          DO j=1, NR2
             CALL rnd_number_bi(rndid2, HARVESTINIPOS2(:,j) &
                  , ng=NGCELL, nl=V_NCELL(0:np-1) )
          END DO

          ! CALLING THE ATTILA MASTER ROUTINE
          CALL ATTILA_INIPOS( &
               HARVESTINIPOS1(:,:), HARVESTINIPOS2(:,:), status)
          IF (status /= 0) &
               CALL error_bi(' ATTILA_INIPOS REPORTED AN ERROR!', substr)

          DEALLOCATE(HARVESTINIPOS1)
          DEALLOCATE(HARVESTINIPOS2)

       ELSE

          !
          ! NOTE: This block contains the new 'single PE'
          !       initialisation, which has the drawback that it is
          !       - for some reason - increadibly slow ...
          !
          CALL rnd_init_bi(rndid1, RND_ATTILA, RND_MP_PIN, ISEED1)
          CALL rnd_init_bi(rndid2, RND_ATTILA, RND_MP_PIN, ISEED2)

          IF (p_pe == p_io) THEN
             ALLOCATE(HARVESTINIPOS1(NGCELL,NR1))
             ALLOCATE(HARVESTINIPOS2(NGCELL,NR2))
          ELSE
             ALLOCATE(HARVESTINIPOS1(0,NR1))
             ALLOCATE(HARVESTINIPOS2(0,NR2))
          END IF

          DO j=1, NR1
             CALL rnd_number_bi(rndid1, HARVESTINIPOS1(:,j))
          END DO
          DO j=1, NR2
             CALL rnd_number_bi(rndid2, HARVESTINIPOS2(:,j))
          END DO

          IF (p_pe == p_io) THEN

             ! SET TEMPORARILY NCELL = NGCELL
             NCELL_SAVE = NCELL
             NCELL      = NGCELL
             ! CREATE TEMPORARY (GLOBAL) SPACE AND SET LOCAL POINTER
             ! ON P_IO TO IT
             CALL ATTILA_REALLOC(1, TMP_APOS, TMP_NPOS, TMP_NPOSM1)
             !
             ! INITIALISE GLOBALLY
             CALL ATTILA_INIPOS( &
                  HARVESTINIPOS1(:,:), HARVESTINIPOS2(:,:), status)
             IF (status /= 0) &
                  CALL error_bi(' ATTILA_INIPOS REPORTED AN ERROR!', substr)

             ! RESET NCELL
             NCELL = NCELL_SAVE
             ! RESET LOCAL POINTERS ON P_IO
             CALL ATTILA_REALLOC(2, TMP_APOS, TMP_NPOS, TMP_NPOSM1 &
                  , LSPAC=LSPAC,LSPACINT=LSPACINT,LSPACINTM1=LSPACINTM1)
          ENDIF

          DEALLOCATE(HARVESTINIPOS1)
          DEALLOCATE(HARVESTINIPOS2)

          ! DISTRIBUTE GLOBAL INITIALISATION FROM P_IO TO ALL PEs
          DO i=1, SIZE(APOS,1)
             IF (p_pe == p_io) THEN
                tmp_g => TMP_APOS(i,:)
             END IF
             tmp_l => APOS(i,:)
             CALL scatter_glix(tmp_g, tmp_l)
          END DO
          DO i=1, SIZE(NPOS, 1)
             IF (p_pe == p_io) THEN
                tmp_g => TMP_NPOS(i,:)
             END IF
             tmp_l => NPOS(i,:)
             CALL scatter_glix(tmp_g, tmp_l)
          END DO
          DO i=1, SIZE(NPOSM1, 1)
             IF (p_pe == p_io) THEN
                tmp_g => TMP_NPOSM1(i,:)
             END IF
             tmp_l => NPOSM1(i,:)
             CALL scatter_glix(tmp_g, tmp_l)
          END DO

          ! CLEAN
          IF (p_pe == p_io) THEN
             ! DEALLOCATE TEMPORARY SPACE ON P_IO
             CALL ATTILA_REALLOC(3, TMP_APOS, TMP_NPOS, TMP_NPOSM1)
          END IF

       END IF l_ini_par
       ! ####################################################################

       CALL rnd_finish_bi(rndid1)
       CALL rnd_finish_bi(rndid2)

!!$    END IF init  ! op_pj_20160803

       ! UPDATE CELL DISTRIBUTION
       CALL attila_update_celldist

    END IF init ! op_pj_20160803

  END SUBROUTINE attila_initialize_positions
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE attila_initialize_positions_t(flag)

    ! - flag=1 at start of attila_global_end:
    !          apos(3) at lstart not correct, if VER_VEL=2
    ! - flag=2 at end of attila_global_end
    ! - called in attila_global_end

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,       ONLY: p_pe, p_bcast, p_nprocs
    ! MESSy
    USE messy_main_timer,        ONLY: lresume, lstart &
                                     , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
    USE messy_main_timer,        ONLY: julian_day

    IMPLICIT NONE
    INTRINSIC :: ALL, NINT

    ! I/O
    INTEGER, INTENT(IN) :: flag

    CHARACTER(LEN=*), PARAMETER :: substr = 'attila_initialize_positions_t'
    INTEGER            :: jk, jc
    REAL(DP)           :: jdf_now   ! julian day (plus fraction) now
    REAl(DP)           :: jdf_trs   ! julian day (plus fraction) of traj. start
    LOGICAL, DIMENSION(:), POINTER :: l_init_finished => NULL()
    LOGICAL, SAVE      :: l_return = .FALSE.

    ! DO NOTHING IF ALL TRAJECTORIES (ON ALL PEs!) HAVE BEEN INITIALIZED
    ! AND ARE BEING INTEGRATED
    IF (l_return) RETURN

    ALLOCATE (l_init_finished(0:p_nprocs-1))
    l_init_finished(p_pe) = ALL(NINT(RTRAJEC_INTEGRATE(:)) >= 0)
    DO jc=0, p_nprocs-1
       CALL p_bcast(l_init_finished(jc), jc)
    END DO
    l_return = ALL(l_init_finished)
    IF (ASSOCIATED(l_init_finished)) THEN
       DEALLOCATE(l_init_finished)
       NULLIFY(l_init_finished)
    END IF

    IF (l_return) RETURN

    flag1: IF (flag == 1) THEN

       first: IF (lstart .OR. lresume) THEN

          ! set positions in SMCL at very first timestep or after restart
          DO jk = 1, NCELL
             IF (RTRAJEC_INTEGRATE(jk) < 0._dp) THEN
                CALL ATTILA_RESETPOS_TRAJEC(jk               &
                     ,90.0_dp - AP(IDX(p_pe,1)-1+jk)%lat     &
                     ,AP(IDX(p_pe,1)-1+jk)%lon               &
                     ,AP(IDX(p_pe,1)-1+jk)%press             )
             ENDIF
          ENDDO
       ENDIF first

    ENDIF flag1

    ! THIS POINT IS ONLY REACHED, IF NOT ALL TRAJECTORIES ARE
    ! BEING INTEGRATED

    ! SET CURRENT DATE
    jdf_now = julian_day(REAL(DAY, DP), MONTH, YEAR) &
         + REAL(HOUR, DP) / 24.0_DP &
         + REAL(MINUTE, DP) / (60.0_DP * 24.0_DP) &
         + REAL(SECOND, DP) / (60.0_DP * 60.0_DP * 24.0_DP)

    cell_loop: DO jk = 1, NCELL

       ! GET START TIME OF TRAJECTORY
       jdf_trs = julian_day(REAL(AP(IDX(p_pe,1)-1+jk)%day, DP) &
            , AP(IDX(p_pe,1)-1+jk)%month &
            , AP(IDX(p_pe,1)-1+jk)%year) &
            + AP(IDX(p_pe,1)-1+jk)%hour / 24.0_DP &
            + AP(IDX(p_pe,1)-1+jk)%minute / (60.0_DP * 24.0_DP)
            ! second = 0

       time: IF (jdf_now >= jdf_trs) THEN

          ! START DATE FOR THIS TRAJECTORY (ALREADY) REACHED
          RTRAJEC_INTEGRATE(jk) = RTRAJEC_INTEGRATE(jk) + 1.0_dp

          IF (LLTINFO) THEN ! op_pj_20160603
             IF (NINT(RTRAJEC_INTEGRATE(jk)) == 0) THEN
                WRITE(*,*) substr,' (p_pe=',p_pe        &
                     , '): START INTEGRATION FOR CELL ' &
                     , IDX(p_pe,1)-1+jk
             ENDIF
          END IF

       ELSE

          ! START DATE FOR THIS TRAJECTORY NOT YET REACHED
          ! -> RESET TO INITIAL POINT
!          WRITE(*,*) substr,' (p_pe=',p_pe        &
!               , '): RESETTING POSITION OF CELL ' &
!               , IDX(p_pe,1)-1+jk

          ! CALL THE ROUTINE DOING THE RESET
          CALL ATTILA_RESETPOS_TRAJEC(jk               &
               ,90.0_dp - AP(IDX(p_pe,1)-1+jk)%lat     &
               ,AP(IDX(p_pe,1)-1+jk)%lon               &
               ,AP(IDX(p_pe,1)-1+jk)%press             )

       END IF time

    END DO cell_loop

    ! BEFORE START OF CALCULATION IN TRAJECTORY MODE,
    ! OUTPUT INIT-POSITIONS INSTEAD OF ZERO
    latz(:)        = 90._dp-LSPAC(2,:,1,1,1)

    ! UPDATE CELL DISTRIBUTION
    CALL attila_update_celldist

  END SUBROUTINE attila_initialize_positions_t
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE attila_integrate

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,       ONLY: p_pe, p_bcast
    USE messy_main_transform_bi, ONLY: trp_gpdc_gpgl
    USE messy_main_timer,        ONLY: lstart

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'attila_integrate'
    ! TI2 (GLOBAL FIELD)
    REAL(dp), DIMENSION(:,:,:), POINTER  :: TI2 => NULL()
    !
    INTEGER :: jk, i, j

    ! ADDITIONAL CONDITIONS: SKIP IN FIRST TIME STEP
    LLCONV = LLCONV_SAVE .AND. (.NOT. lstart)
    LLCAT  = LLCAT_SAVE  .AND. (.NOT. lstart)

    IF (LLCAT) CALL trp_gpdc_gpgl(1, cati, TI2)

    DO i=1, NRPROC
       proc_active: IF (LPROC(i)) THEN
! qqq TEMPORAL SPECIAL CASE FOR CONV
          IF (i == RANDI_CONV) THEN
             ! no jump ahead, same sequence of length NGCELL on all PEs
             DO j=1, RAND_MULT(i)
                CALL rnd_number_bi(rndid(i), HARVEST(i)%ptr(:,j))
             END DO
          ELSE
             ! jump ahead
             DO j=1, RAND_MULT(i)
                ! NOTE: Parameters ng AND nl are obsolete, if
                !       I_RANDOM_PARALLEL = RND_MP_PIJ
                CALL rnd_number_bi(rndid(i), HARVEST(i)%ptr(:,j) &
                     , ng=NGCELL, nl=V_NCELL(0:np-1) )
             END DO
          END IF
          IF (L_RANDOM_TEST .AND. (i == RANDI_TURB)) THEN
             test_rnd(:) = HARVEST(i)%ptr(:,1)
          END IF
       ELSE
          ! NOTE:
          !   - 2*<random>-1 e [-1,1] => 2*0.5-1=0 for LLTMOCA
          !   0 for LLTBLTURB, LCAT, LCONV
          HARVEST(i)%ptr(:,:) = RDEF(i)  ! default
       END IF proc_active
       !
    END DO

    ! NOTE:
    !  - RANDI_TURB                            1 -> PBL-TURBULENCE
    !  - RANDI_MOCA                            4 -> MC-DIFFUSION
    !  - RANDI_CONV        -> ATTILA_CONVTRAJ  2 -> CONVECTION
    !  - RANDI_CAT   & TI2                     3 -> CLEAR AIR TURBULENCE
    !
    CALL ATTILA_DRIVE( &
           HARVEST(RANDI_TURB)%ptr(:, 1) &
         , HARVEST(RANDI_MOCA)%ptr(:, :) &
         , HARVEST(RANDI_CAT )%ptr(:, 1) &
         , TI2, LSTART)

    ! UPDATE CELL DISTRIBUTION
    CALL attila_update_celldist

    IF (LTRAJEC) THEN
       ! TO ENSURE THAT THE PRESSURE HEIGHT OF ALL AIR PARCELS THAT HAVE
       ! NOT BE INITIALIZED YET REMAINS UNCHANGED IN THE OUTPUT FILE, RESET THE
       ! CURRENT VALUE AND SAVE IN THE CHANNEL THE INITIAL PRESSURE HEIGHT
       !
       ! THIS MUST BE DONE HERE, SINCE ATTILA_COUNT_CELLS ABOVE ADJUSTS
       ! THE PRESSURE ACCORDING TO THE SURFACE PRESSURE IN EVERY TIME STEP
       DO jk=1, NCELL
          IF (RTRAJEC_INTEGRATE(jk) >= 0) CYCLE
          pressheight(jk)=AP(IDX(p_pe,1)-1+jk)%press
       END DO
    END IF

    ! UPDATE GEOGRAPHICAL LATITUDE ET. AL
    latz(:) = 90.0_dp - LSPAC(2,:,1,1,1)

    ! COPY MASS OF AIR PARCEL TO CHANNEL
    ptramcell=AMCELL

    ! RE-SET CONDITIONS
    LLCONV = LLCONV_SAVE
    LLCAT  = LLCAT_SAVE

    ! CLEAN UP
    IF (ASSOCIATED(TI2)) DEALLOCATE(TI2)

  END SUBROUTINE attila_integrate
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE attila_update_celldist

    ! ECHAM5/MESSy
    USE messy_main_transform_bi, ONLY: trp_gpdc_gpgl, M_SUM

! ob_sb_20140516+
    USE messy_main_timer,           ONLY: lstart !!$, lresume
! ob_sb_20140516-

    IMPLICIT NONE

    CALL attila_count_cells(NCB_P, NCBL_P &
         , lstart .AND. .NOT. L_INI_FINI)

    ! SUM NUMBER OF CELLS OF ALL CPUs
    CALL trp_gpdc_gpgl (-1, SGNCB,  NCB_P,  M_SUM)
    CALL trp_gpdc_gpgl (-1, SGNCBL, NCBL_P, M_SUM)

    ! MAKE NUMBER OF CELLS PER GRID-BOX AVAILABLE ON ALL CPUs
    ! (needed for MASS CONSERVING TRANSFORMATION LG <-> GP)
    CALL trp_gpdc_gpgl(1, SGNCB,  GNCB)
    CALL trp_gpdc_gpgl(1, SGNCBL, GNCBL)

  END SUBROUTINE attila_update_celldist
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE attila_convect

    ! ECHAM5/MESSy
    USE messy_main_transform_bi,    ONLY: trp_gpdc_gpgl, M_SUM
    USE messy_main_mpi_bi,          ONLY: p_pe
    USE messy_main_timer,           ONLY: lstart

    IMPLICIT NONE

    ! I/O
    INTEGER :: status

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'attila_convect'
    REAL(DP), DIMENSION(NCELL)  :: kstart
    LOGICAL :: L_EXIT

    IF (lstart) RETURN

    ! CONVECTION TYPE, UPDRAFT AND DOWNDRAFT MASS FLUX [KG/M^2/S]
    CALL trp_gpdc_gpgl(1, uflx_conv, umassfl)
    CALL trp_gpdc_gpgl(1, dflx_conv, dmassfl)
    CALL trp_gpdc_gpgl(1, type_conv, typeconv)

    CALL ATTILA_CONVGPC(umassfl,dmassfl,typeconv,status)
    IF (status /= 0) CALL error_bi('error in ATTILA_CONVGPC', substr)

    ! mz_pj_20080131+
    ! NOTE: For the calculation of convective updraft/downdraft/subsidence
    !       with the usual approach of independent random number vectors
    !       2*2*NLEV+2 random number vectors of length NGCELL are required.
    !       Most entries, however, are not used and the required memory
    !       corresponds to CPGBAVE*(4*NLEV+2) 3D grid-point fields (!!!).
    !       Therefore another approach with only one random number vector
    !       is followed here: Whenever a random number is required, we
    !       shift the random number position by one. In order to make the
    !       result decomposition independent, all PEs have to get the full
    !       random number vector such that the CELL number can
    !       unambiguously be associated with the random number position.
    !       Recycling the random numbers, however, can in principle cause
    !       undesired correlations ... This has to be checked!
    ! mz_pj_20080131-
!qqq THIS SHOULD BE REVISED: USE STANDARD SEQUENCE OF LENGTH NCELL,
    ! COUNT USED NUMBERS AND HARVEST NEW BLOCK IN CASE ALL NUMBERS HAVE
    ! BEEN CONSUMED. IMPLEMENTATION NEEDS TO BE DECOMP-INDEPENDENT!!!

    CALL attila_save_pos(kstart, 3)

    ! ------------
    ! up/downdraft
    ! ------------
    CALL attila_convtraj(1, IDX(p_pe,1)-1, HARVEST(RANDI_CONV)%ptr(:, 1) &
         , POS_BEG_UP, POS_END_UP, POS_END_DO, POS_END_SUB, L_EXIT)

    ! COUNT CELLS ON EACH PE (after up/downdraft)
    CALL attila_count_cells(NCB_P, NCBL_P)

    ! SUM NUMBER OF CELLS OF ALL PEs (after up/downdraft)
    CALL trp_gpdc_gpgl (-1, SGNCB_CUD,  NCB_P,  M_SUM)

    ! MAKE NUMBER OF CELLS PER GRID-BOX AVAILABLE ON ALL PEs
    ! NOTE: THIS IS THE NUMBER OF CELLS TO BE MOVED BY SUBSIDENCE
    !       IN ORDER TO RE-ESTABLISH THE CELL DISTRIBUTION BEFORE CONVECTION
    CALL trp_gpdc_gpgl(1, SGNCB_CUD,  GNCB_MOVE)

    ! ----------
    ! subsidence
    ! ----------
    NCONV_ITER = 0.0_dp

    DO
       IF (L_EXIT) EXIT

       CALL attila_convtraj(2, IDX(p_pe,1)-1, HARVEST(RANDI_CONV)%ptr(:, 1) &
            , POS_BEG_UP, POS_END_UP, POS_END_DO, POS_END_SUB, L_EXIT)

       ! COUNT CELLS ON EACH PE (after subsidence)
       CALL attila_count_cells(NCB_P, NCBL_P)

       ! SUM NUMBER OF CELLS OF ALL PEs (after subsidence)
       CALL trp_gpdc_gpgl (-1, SGNCB_CSUB,  NCB_P,  M_SUM)

       ! MAKE NUMBER OF CELLS PER GRID-BOX AVAILABLE ON ALL PEs
       ! NOTE: THE SUBSIDENCE IS TO RE-ESTABLISH THE CELL DISTRIBUTION
       !       BEFORE CONVECTION; THUS DURING THE ITERATION, SGNCB_CSUB
       !       MUST CONVERGE TO SGNCB ...
       CALL trp_gpdc_gpgl(1, SGNCB_CSUB,  GNCB_MOVE)

       NCONV_ITER = NCONV_ITER + 1.0_dp

    END DO

    CALL attila_save_pos(CONVDELK, 3)
    CONVDELK(:) = CONVDELK(:)-kstart(:)

  END SUBROUTINE attila_convect
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE attila_calc_sidot

 ! ECHAM5/MESSy
  USE messy_main_data_bi,     ONLY: pressi_3d,     &
                                    aps,            &
                                    vervel_scb, alpste_scb
  USE messy_main_grid_def_mem_bi, ONLY: nproma, ngpblks ,      &
                                        npromz , nlev

  IMPLICIT NONE

  ! LOCAL
  INTEGER :: jl, jk, jg, zproma

  sidot(:,:,:) = 0._dp

  DO jg = 1, ngpblks

      if ( jg == ngpblks) then
         zproma = npromz
      else
         zproma = nproma
      endif

      do jl =1, zproma

          do jk=1,nlev
          sidot(jl,jk,jg) = vervel_scb(jl,jk,jg)/aps(jl,jg) - &
               pressi_3d(jl,jk,jg)  &
               /aps(jl,jg) * alpste_scb(jl,jg)
          enddo
      enddo
  enddo

  END SUBROUTINE attila_calc_sidot
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE attila_calc_xidot

  ! ECHAM5/MESSy
  USE messy_main_data_bi,         ONLY: pressi_3d, press_3d, tpot_3d, &
                                        etadot_3d,                    &
                                        aps, alpste_scb, vervel_scb
  USE messy_main_grid_def_mem_bi, ONLY: apzero
  USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev, ngpblks,        &
                                        npromz, nlevp1, nlevm1,       &
                                        vct

  USE messy_main_constants_mem,  ONLY: pi

  IMPLICIT NONE
  INTRINSIC :: SIN, COS

  ! LOCAL
  REAL(dp) :: ref_func(nproma, nlevp1,ngpblks) !op_sb_20160510
  REAL(dp) :: dref_func(nproma, nlevp1)
  REAL(dp) :: sigma(nproma, nlevp1)
  REAL(dp) :: spressi(nproma, nlevp1)
  REAL(dp) :: dpressi(nproma, nlev)
  REAL(dp) :: pit2

  INTEGER :: jl, jk, jg, zproma
  INTEGER :: marker(nproma,ngpblks)

  ALLOCATE(sigma_ref(nproma,ngpblks))
  ALLOCATE(i_ref(nproma,ngpblks))

  xidot(:,:,:)    = 0.0_dp
  dthetadt(:,:,:) = 0.0_dp
  spressi(:,:)    = 0._dp
  dpressi(:,:)    = 0._dp
  sigma(:,:)      = 0._dp
  ref_func(:,:,:) = 0._dp
  dref_func(:,:)  = 0._dp
  sigma_ref(:,:)  = 0._dp
  i_ref(:,:)      = 0._dp

  pit2 = pi / 2._dp

  DO jg = 1, ngpblks

      if ( jg == ngpblks) then
         zproma = npromz
      else
         zproma = nproma
      endif

      do jl =1, zproma

         sigma(jl, :)    = pressi_3d(jl,:,jg)  / aps(jl,jg)
         sigma(jl,nlevp1) = MIN( sigma(jl,nlevp1), 1._dp )

         if (press_ref < 0.0_dp) then
           sigma_ref(jl,jg)  = tp(jl,jg) / apzero         !aps(jl,jg)
         else
           sigma_ref(jl,jg)  = press_ref / apzero         !aps(jl,jg)
         endif

         marker(jl,jg)=0
         do jk = 1, nlevp1
            if (pressi_3d(jl,jk,jg) >= (sigma_ref(jl,jg)*apzero)) then
               ! troposphere (below ca. 100 hPa)
               ref_func(jl,jk,jg) = sin(pit2 * &
                    (1._dp-pressi_3d(jl,jk,jg)/aps(jl,jg)) &
                    / (1._dp-sigma_ref(jl,jg)))
               if (marker(jl,jg) == 0) then
                   i_ref(jl,jg) = real(jk-1)
                   marker(jl,jg) = 1
               endif

                dref_func(jl,jk) = pit2*cos(pit2 * (1._dp-sigma(jl,jk)) &
                                     / (1._dp-sigma_ref(jl,jg))) *   &
                                   (-sigmadot(jl,jk,jg)/ &
                                   (1._dp-sigma_ref(jl,jg)))
               ! op_pj_20150310+: avoid numerical problems;
               ref_func(jl, jk, jg) = MAX(ref_func(jl, jk, jg), 0.0_dp)
               ! op_pj_20150310-
            else
               ! stratosphere+above 100 hPa)
               ref_func(jl,jk,jg)  = 1._dp
               dref_func(jl,jk)    = 0._dp
            endif
         enddo

         do jk=1,nlev
            ! calculate dtdt for diabatic processes only
            dpressi(jl,jk)     = pressi_3d(jl,jk+1,jg)-pressi_3d(jl,jk,jg)

            ! diabatic heating rates Q*theta/T
            dthetadt(jl,jk,jg) = tempdot(jl,jk,jg) *     &
                                (100000._dp/press_3d(jl,jk,jg))**kappa

         enddo

         do jk=1,nlevm1
            spressi(jl,jk) =  dpressi(jl,jk+1) + dpressi(jl,jk)
         enddo

         xidot(jl,1,jg) = 0._dp
         xidot(jl,nlevp1,jg) = 0._dp

         do jk=2,nlev

            diab(jl,jk,jg) =  ref_func(jl,jk,jg)*(dthetadt(jl,jk-1,jg)*   &
                 dpressi(jl,jk)/spressi(jl,jk-1)+dthetadt(jl,jk,jg)* &
                 dpressi(jl,jk-1) /spressi(jl,jk-1))

            kine(jl,jk,jg) =  dref_func(jl,jk) *    &
                 (tpot_3d(jl,jk-1,jg)* dpressi(jl,jk)/spressi(jl,jk-1) &
                 +tpot_3d(jl,jk,jg) * dpressi(jl,jk-1)/spressi(jl,jk-1))

            xidot(jl,jk,jg)=diab(jl,jk,jg) + kine(jl,jk,jg)

            diab(jl,jk,jg) = ref_func(jl,jk,jg)
            kine(jl,jk,jg) = dref_func(jl,jk)

         enddo
      enddo
   END DO

  END SUBROUTINE attila_calc_xidot
!---------------------------------------------------------------------------

! ===========================================================================
END MODULE messy_attila_e5
! ===========================================================================
