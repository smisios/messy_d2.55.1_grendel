#include "messy_main_ppd_bi.inc"

MODULE  messy_ddep_si

  ! BML/MESSy
  USE messy_main_tracer_mem_bi, ONLY: ntrac_gp, ntrac_lg, ti_gp, ti_lg &
                                    , GPTRSTR,  LGTRSTR
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                    , info_bi, warning_bi, error_bi
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle,       &
                                      mtend_get_start_l,      &
                                      mtend_add_l,            &
                                      mtend_register,         &    
                                      mtend_id_tracer
#endif

  ! MESSy
  USE messy_main_constants_mem, ONLY: g, M_air, N_A
  USE messy_main_channel,       ONLY: STRLEN_OBJECT, STRLEN_CHANNEL     &
                                    , t_chaobj_cpl
  USE messy_main_tracer,        ONLY: ON, I_drydep, R_molarmass, R_pss  &
                                    , R_dryreac_sf, I_Aerosol_mode      &
                                    , S_aerosol_model, I_Aerosol_method &
                                    , get_tracer                        &
                                    , STRLEN_FNAME
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY, PTR_2D_ARRAY &
                                    , PTR_4D_ARRAY, PTR_5D_ARRAY
  ! SUBMODEL DDEP
  USE messy_ddep

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE
  SAVE

  REAL(DP), PARAMETER :: conv_mr = N_A / (1.E-03_dp * M_air)
  REAL(DP), PARAMETER :: conv_nu = 1._dp / (1.E-03_dp * M_air)

  TYPE LOGICAL_1D_ARRAY
    LOGICAL,  DIMENSION(:), POINTER :: LPTR => NULL()
  END TYPE LOGICAL_1D_ARRAY

  ! NEW TYPES FOR DEPOSITION WITH PRESCRIBED DEPOSITION VELOCITIES

  ! NMAX_PREDEPVEL: MAX. NUMBER OF PREDEFINED DEPOSITION VELOCITY FIELDS
  INTEGER, PARAMETER :: NMAX_PREDEPVEL = 250
  ! ACTUAL NUMBER OF VELOCITY FIELDS / TRACERS
  INTEGER            :: N_predepvel 

  TYPE T_IO_PREDEPVEL 
     ! channel of predefined deposition velocity field (usually import_grid)
     CHARACTER(LEN=STRLEN_CHANNEL) :: CHANNEL = ''     
     ! name of predefined deposition velocity field
     CHARACTER(LEN=STRLEN_OBJECT)  :: OBJECT = ''     
     ! string defining the tracer being affected by the predefined 
     ! deposition velocity field
     CHARACTER(LEN=STRLEN_FNAME)   :: TRNAME = ''     
     ! dry-deposition flux diagnostics 
     LOGICAL                       :: TRLOGICAL = .FALSE. 
  END TYPE T_IO_PREDEPVEL 

  TYPE(T_IO_PREDEPVEL), DIMENSION(NMAX_PREDEPVEL) :: import_predepvel

  TYPE T_PREDEPVEL  
     TYPE(T_IO_PREDEPVEL)                :: import_predepvel
     ! ptr to predefined deposition velocity field 
     REAL(DP), DIMENSION(:,:), POINTER   :: field_predepvel  => NULL()
     ! diagnostic ptrs
     REAL(DP), DIMENSION(:,:), POINTER   :: diagflux  => NULL()
     REAL(DP), DIMENSION(:,:), POINTER   :: diagsum  => NULL()
     ! name of tracer subject to predefined dry deposition)
     CHARACTER(LEN=STRLEN_MEDIUM)        :: NAME = '' 
     ! subname of tracer subject to predefined dry deposition
     CHARACTER(LEN=STRLEN_MEDIUM)        :: SUBNAME = ''
     ! ID of tracer subject to predefined dry deposition
     INTEGER                             :: trid = 0
  END TYPE T_PREDEPVEL

  TYPE(T_PREDEPVEL), DIMENSION(NMAX_PREDEPVEL) :: X_predepvel    

  ! CPL-NAMELIST PARAMETERS
  LOGICAL :: L_LG = .FALSE.  ! dry deposition for Lagrangian tracers
  LOGICAL :: L_GP = .FALSE.  ! dry deposition for Gridpoint tracers
  INTEGER :: i_lg_method = 1 ! dry deposition method (LG)
  ! switch for writing change directly in tendency 
  LOGICAL, PUBLIC :: l_tendency = .TRUE.
  ! switch for diagnostic output of deposition fluxes
  LOGICAL, PUBLIC :: l_diagflux = .FALSE.
  ! extra output of dry deposition fluxes for special tracers required
  LOGICAL, PUBLIC            :: l_fluxout = .FALSE. !qqq _gp, _lg
  CHARACTER(LEN=1024), PUBLIC:: outflux =''
  TYPE(t_chaobj_cpl)         :: imp_lai
  TYPE(t_chaobj_cpl)         :: imp_hc
  TYPE(t_chaobj_cpl)         :: imp_drag
  TYPE(t_chaobj_cpl)         :: imp_soilpH

  TYPE(t_chaobj_cpl)         :: rainrate_ls
  TYPE(t_chaobj_cpl)         :: rainrate_cv

  ! DATASET NUMBERING
  INTEGER, PARAMETER :: I_GP = 1           ! GP is always DSET 1
  INTEGER, PARAMETER :: I_LG = 2           ! LG is always DSET 2
  CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: setname = (/GPTRSTR, LGTRSTR/)

  !  FOR AEROSOL DEPOSITION
  INTEGER, PARAMETER :: NMAX_AEROMODELS = 10

  ! DEFINE DATA STRUCT FOR GP AND LG TRACER
  TYPE T_DRYDEP_SET
     INTEGER  :: reprid = 0                   ! representation ID OF SET 
     LOGICAL  :: l_vdbl   = .FALSE.           ! big leave gas phase 
     LOGICAL  :: l_vdaer  = .FALSE.           ! calc. aerosol dry deposition
     LOGICAL  :: l_airsea_coupling = .FALSE.
     LOGICAL, DIMENSION(:), POINTER  :: l_cpl_airsea => NULL()
     INTEGER  :: nodoc = 0             ! number of diagnostic output channels
     INTEGER  :: ntrac = 0
     INTEGER, POINTER, DIMENSION(:) :: idt_vdbl => NULL()
     INTEGER, POINTER, DIMENSION(:) :: idt_diag => NULL() 
     INTEGER                        :: idt_h2so4
     TYPE (PTR_2D_ARRAY), DIMENSION(:), POINTER :: drydepflux         => NULL()
     TYPE (PTR_2D_ARRAY), DIMENSION(:), POINTER :: drydepfluxsum      => NULL()
     TYPE (PTR_2D_ARRAY), DIMENSION(:), POINTER :: vd_bigl            => NULL()
     TYPE (PTR_2D_ARRAY), DIMENSION(:), POINTER :: drydepflux_lggp    => NULL()
     TYPE (PTR_2D_ARRAY), DIMENSION(:), POINTER :: drydepfluxsum_lggp => NULL()
     REAL(dp), DIMENSION(:,:,:), POINTER        :: zvdrydep    => NULL()
     REAL(dp), DIMENSION(:,:),   POINTER        :: zvdrydep_lg => NULL()
     REAL(dp), DIMENSION(:,:,:), POINTER        :: vd_eff      => NULL()
     REAL(dp), DIMENSION(:,:),   POINTER        :: vd_eff_lg   => NULL()
     ! FIELDS FOR AEROSOL DEPOSITION
     CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(NMAX_AEROMODELS) :: aermodname 
     INTEGER,                      DIMENSION(NMAX_AEROMODELS) :: aermodmethod
     INTEGER                                                  :: aeromodnum = 0
     ! AEROSOL MODEL (ID) OF TRACER
     INTEGER, DIMENSION(:), POINTER             :: traermod => NULL()
     TYPE (PTR_2D_ARRAY), DIMENSION(:), POINTER :: radii_lg   => NULL()
     TYPE (PTR_2D_ARRAY), DIMENSION(:), POINTER :: densaer_lg => NULL()
     TYPE (PTR_4D_ARRAY), DIMENSION(:), POINTER :: radii   => NULL()
     TYPE (PTR_4D_ARRAY), DIMENSION(:), POINTER :: densaer => NULL()
     TYPE (PTR_1D_ARRAY), DIMENSION(:), POINTER :: sigma   => NULL()
     INTEGER,             DIMENSION(:), POINTER :: nmod => NULL()

     TYPE (LOGICAL_1D_ARRAY), DIMENSION(:), POINTER :: l_anymode => NULL()
     TYPE (LOGICAL_1D_ARRAY), DIMENSION(:), POINTER :: l_anymass => NULL()
     TYPE (LOGICAL_1D_ARRAY), DIMENSION(:), POINTER :: l_anynum  => NULL()
     !
     ! deposition velocity for mass distribution
     TYPE (PTR_5D_ARRAY), DIMENSION(:), POINTER :: vdrydep_mass => NULL() 
     ! deposition velocity for number distribution
     TYPE (PTR_5D_ARRAY), DIMENSION(:), POINTER :: vdrydep_num  => NULL()
     !
     REAL(dp), DIMENSION(:), POINTER :: diff => NULL()
     REAL(dp), DIMENSION(:), POINTER :: diffrb => NULL() 
     REAL(dp), DIMENSION(:), POINTER :: rmes => NULL()
     REAL(dp), DIMENSION(:), POINTER :: rcut => NULL()
     ! ju_te_20180625+
     TYPE (PTR_2D_ARRAY), DIMENSION(:), POINTER :: rws_2d  => NULL()
     TYPE (PTR_2D_ARRAY), DIMENSION(:), POINTER :: rcut_2d => NULL()
     ! ju_te_20180625-
     REAL(dp), DIMENSION(:), POINTER :: rsoil => NULL() 
     REAL(dp), DIMENSION(:), POINTER :: rws => NULL()
     REAL(dp), DIMENSION(:), POINTER :: rwater => NULL()
     REAL(dp), DIMENSION(:), POINTER :: rsnow => NULL()
  END TYPE T_DRYDEP_SET

  TYPE(T_DRYDEP_SET), DIMENSION(:), POINTER :: DSET
  !
  ! INPUT DATEN
  REAL(dp), DIMENSION(:,:),   POINTER :: lai  => NULL()     ! LAI
  REAL(dp), DIMENSION(:,:),   POINTER :: hc   => NULL()     ! canopy height
  REAL(dp), DIMENSION(:,:),   POINTER :: drag => NULL()     ! drag
  REAL(dp), DIMENSION(:,:),   POINTER :: z0m  => NULL()     ! surface roughness
  REAL(dp), DIMENSION(:,:,:), POINTER :: soilph => NULL()   ! soil pH classes
  !
  ! large-scale rain rate ('kg m-2 s-1')
  REAL(dp), DIMENSION(:,:),   POINTER :: rrain_ls  => NULL() 
  ! convective rain rate  ('kg m-2 s-1')
  REAL(dp), DIMENSION(:,:),   POINTER :: rrain_cv  => NULL() 

  ! COUPLING TO ATTILA SUBMODEL
  ! AIR MASS OF CELL
  REAL(dp), POINTER               :: MASS_PARCEL => NULL()
  ! PRESS_PARCEL(NCELLS) = pressure coordinate of the parcel
  REAL(dp), POINTER, DIMENSION(:) :: PRESS_PARCEL => NULL() 
 
#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  ! SUBROUTINES
  PUBLIC :: ddep_initialize
  PUBLIC :: ddep_init_memory
  PUBLIC :: ddep_init_coupling
  PUBLIC :: ddep_global_start
  PUBLIC :: ddep_vdiff
  PUBLIC :: ddep_global_end
  PUBLIC :: ddep_free_memory

  !PRIVATE :: calc_vd_drydep
  !PRIVATE :: drydep_read_nml_cpl
  !PRIVATE :: ddep_predef_vdiff

CONTAINS

  !==========================================================================

  SUBROUTINE ddep_initialize

    ! DDEP MODULE ROUTINE (SUBMODEL INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    ! INITIALIZATION OF DRYDEP SPECIFIC EVENTS FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002       (DRYDEP)
    ! Modified, Laurens Ganzeveld, MPICH, 31-10-2002 (DRYDEP)
    ! Modified, Astrid Kerkweg, MPICH, Mar 2004      (DRYDEP)
    ! 
    ! Astrid Kerkweg, UNI-MZ, Apr 2010 new submodel DDEP without RGT

    ! BML/MESSy
    USE messy_main_mpi_bi,         ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,          ONLY: find_next_free_unit &
                                       , strcrack
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
  
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='ddep_initialize'
    INTEGER                     :: status
    INTEGER                     :: iou    ! I/O unit

    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: &
         strcrack_tmp => NULL()
    INTEGER                     :: i, j, strcrack_size

    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL DDEP CORE ROUTINE:
       CALL ddep_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    ! BROADCAST RESULTS
    CALL p_bcast(l_whitecap, p_io)
    CALL p_bcast(l_rh, p_io)
    CALL p_bcast(l_ganzeori, p_io) !ju_te_20180625

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL DDEP CORE ROUTINE:
       CALL ddep_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    ! BROADCAST RESULTS
    CALL p_bcast(l_lg, p_io)
    CALL p_bcast(l_gp, p_io)

    CALL p_bcast(i_lg_method, p_io)
    CALL p_bcast(l_tendency,  p_io)
    CALL p_bcast(l_diagflux,  p_io)
    CALL p_bcast(outflux,     p_io)

    CALL p_bcast(imp_lai%CHA,    p_io)
    CALL p_bcast(imp_lai%OBJ,    p_io)
    CALL p_bcast(imp_drag%CHA,   p_io)
    CALL p_bcast(imp_drag%OBJ,   p_io)
    CALL p_bcast(imp_hc%CHA,     p_io)
    CALL p_bcast(imp_hc%OBJ,     p_io)
    CALL p_bcast(imp_soilpH%CHA, p_io)
    CALL p_bcast(imp_soilpH%OBJ, p_io)

    CALL p_bcast(rainrate_ls%CHA, p_io)
    CALL p_bcast(rainrate_ls%OBJ, p_io)
    CALL p_bcast(rainrate_cv%CHA, p_io)
    CALL p_bcast(rainrate_cv%OBJ, p_io)


    ! predefined dry deposition: import into X_predepvel and broadcast
    IF (p_parallel_io) THEN
       ! copy data and parse trname 
       N_predepvel = 1 ! determine how many surface sinks
       DO i=1, NMAX_PREDEPVEL
          IF (TRIM(import_predepvel(i)%CHANNEL) == '') CYCLE
          X_predepvel(N_predepvel)%import_predepvel%CHANNEL &
               = import_predepvel(i)%CHANNEL
          WRITE(*,*) 'PREDEFINED DEPOSITION VELOCITY ''' &
               , TRIM(X_predepvel(N_predepvel)%import_predepvel%CHANNEL) &
               ,''' ...'
          X_predepvel(N_predepvel)%import_predepvel%OBJECT    = &
               import_predepvel(i)%OBJECT
          X_predepvel(N_predepvel)%import_predepvel%TRNAME    = &
               import_predepvel(i)%TRNAME
          X_predepvel(N_predepvel)%import_predepvel%TRLOGICAL = &
               import_predepvel(i)%TRLOGICAL
          CALL strcrack( X_predepvel(N_predepvel)%import_predepvel%TRNAME, '_' &
               , strcrack_tmp, strcrack_size)
          SELECT CASE (strcrack_size)
          CASE (1)
             X_predepvel(N_predepvel)%NAME = strcrack_tmp(1)
             X_predepvel(N_predepvel)%SUBNAME = ''
          CASE DEFAULT 
             X_predepvel(N_predepvel)%NAME = strcrack_tmp(1)
             X_predepvel(N_predepvel)%SUBNAME = TRIM(strcrack_tmp(2))
             DO j=3, strcrack_size
                X_predepvel(N_predepvel)%SUBNAME = &
                     TRIM(X_predepvel(N_predepvel)%SUBNAME)//'_'//&
                     &TRIM(strcrack_tmp(j))
             END DO
          END SELECT
          N_predepvel = N_predepvel + 1
       END DO
       N_predepvel = N_predepvel - 1
    END IF
    CALL p_bcast(N_predepvel, p_io)
    DO i=1, N_predepvel !loop number of surface sinks
       CALL p_bcast(X_predepvel(i)%import_predepvel%CHANNEL, p_io)
       CALL p_bcast(X_predepvel(i)%import_predepvel%OBJECT,  p_io)
       CALL p_bcast(X_predepvel(i)%import_predepvel%TRNAME,  p_io)
       CALL p_bcast(X_predepvel(i)%import_predepvel%TRLOGICAL,  p_io)
       CALL p_bcast(X_predepvel(i)%NAME,  p_io)
       CALL p_bcast(X_predepvel(i)%SUBNAME, p_io)
    END DO
    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',N_predepvel,' CASES WITH PREDEFINED DRY DEPOSITION INITIALIZED !'
    END IF

  END SUBROUTINE ddep_initialize

! ========================================================================

  SUBROUTINE ddep_init_memory

    ! BML/MESSy
    USE messy_main_mpi_bi,            ONLY: p_parallel_io
    USE messy_main_channel_error_bi,  ONLY: channel_halt
    USE messy_main_channel_bi,        ONLY: GP_2D_HORIZONTAL, LG_ATTILA 
    USE messy_main_tracer_mem_bi,     ONLY: NCELL
    USE messy_main_grid_def_mem_bi,   ONLY: nproma, ngpblks
#ifdef MESSYTENDENCY
    USE messy_main_tracer_tools_bi,   ONLY: tracer_halt
#endif

    ! MESSy
    USE messy_main_tracer,        ONLY: AEROSOL, AIR, t_trinfo_tp &
                                          , NUMBERDENSITY, AMOUNTFRACTION
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_attribute, get_channel_info
    USE messy_main_tools,         ONLY: strcrack, int2str

    IMPLICIT NONE

    INTRINSIC TRIM

    ! LOCAL
    CHARACTER(len=*), PARAMETER  :: substr = 'ddep_init_memory'
    INTEGER                      :: jt, naermo, ii, i
    ! CHANNEL MANAGEMENT
    INTEGER                      :: status
    LOGICAL                      :: lnew
    INTEGER                      :: nvdbl
    INTEGER                      :: na

    TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti => NULL()
    CHARACTER(LEN=26), POINTER   :: strname(:)  => NULL()
    CHARACTER(LEN=2)             :: teststr
    INTEGER                      :: offset, ic
    CHARACTER(LEN=STRLEN_MEDIUM) :: unit
#ifdef MESSYTENDENCY
    INTEGER                      :: idt
#endif

    CALL start_message_bi(modstr, 'MEMORY INITIALIZATION', substr)

    !------------------------------------------------------------------------
    !   GENERAL PART USED BY LAGRANGE AND GRID POINT PART
    !------------------------------------------------------------------------
#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif
    
    ! (1) TEST AND ALLOCATE NUMBER OF DRYDEP TRACER SETS

    IF (L_GP .AND. ntrac_gp == 0) THEN
       CALL warning_bi('L_GP is .TRUE., but no GP-Tracer available', substr)
       CALL warning_bi('SETTING L_GP to .FALSE.', substr)
       L_GP = .FALSE.
    ENDIF
#ifdef ECHAM5
    IF (L_LG .AND. ntrac_lg == 0) THEN
!!#D attila +
       CALL warning_bi('L_LG is .TRUE., but no LG-Tracer available', substr)
       CALL warning_bi('SETTING L_LG to .FALSE.', substr)
!!#D attila -
       L_LG = .FALSE.
    ENDIF
#endif

    IF (L_GP .AND. L_LG) THEN
       ALLOCATE(DSET(2))
    ELSE IF (L_GP) THEN
       ALLOCATE(DSET(I_GP:I_GP))
    ELSE IF (L_LG) THEN
       ALLOCATE(DSET(I_LG:I_LG))
    ELSE
       CALL warning_bi('NO TRACER AVAILABLE NEITHER GP NOR LG',substr)
       CALL error_bi('PLEASE SWITCH OF DRYDEP or implement new tracer set'&
            ,substr)
       RETURN
    ENDIF

    ! (2) INTERNAL UPDATE OF SWITCHES

    IF (L_GP) THEN
       DO jt = 1, ntrac_gp
          IF (ti_gp(jt)%tp%meta%cask_i(I_drydep) == ON) &
                                                DSET(I_GP)%l_vdbl = .TRUE.
          IF (ti_gp(jt)%tp%ident%medium == AEROSOL .AND. &
               ti_gp(jt)%tp%meta%cask_i(I_drydep) == ON ) &
                                               DSET(I_GP)%l_vdaer = .TRUE.
       END DO
       IF (.NOT. DSET(I_GP)%l_vdbl .AND. DSET(I_GP)%l_vdaer) &
            CALL info_bi( &
            ' NO GASPHASE TRACER IN GRID POINT SPACE with ndrydep==ON, '//&
            &'but Aerosol', substr )
       
       IF (DSET(I_GP)%l_vdbl) &
            CALL info_bi('big leaf dry deposition switched on for GP', '  ')
       IF (DSET(I_GP)%l_vdaer) &
            CALL info_bi ('aerosol dry deposition switched on for GP', '  ')

       IF (.NOT. DSET(I_GP)%l_vdbl .AND. .NOT. DSET(I_GP)%l_vdaer) THEN
          CALL info_bi('NO TRACER IN GRID POINT SPACE with ndrydep=ON',substr)
          CALL info_bi('SWITCHING OF  L_GP',substr)
          L_GP = .FALSE.
       ENDIF
    ENDIF

!!#D attila +
#ifdef ECHAM5
    IF (L_LG) THEN
       DO jt = 1, ntrac_lg
          IF (ti_lg(jt)%tp%meta%cask_i(I_drydep) == ON) &
                                                DSET(I_LG)%l_vdbl = .TRUE.
          IF (ti_lg(jt)%tp%ident%medium == AEROSOL .AND. &
               ti_lg(jt)%tp%meta%cask_i(I_drydep) == ON ) &
                                               DSET(I_LG)%l_vdaer = .TRUE.
       END DO
       IF (.NOT. DSET(I_LG)%l_vdbl .AND. DSET(I_LG)%l_vdaer) &
            CALL info_bi( &
            ' NO GASPHASE TRACER IN LAGRANGIAN SPACE with ndrydep==ON, '//&
            &'but Aerosol',substr )

       IF (DSET(I_LG)%l_vdbl) &
            CALL info_bi('big leaf dry deposition switched on for LG', '  ')
       IF (DSET(I_LG)%l_vdaer) &
            CALL info_bi('aerosol dry deposition switched on for LG', '  ')
       IF (.NOT. DSET(I_LG)%l_vdbl .AND. .NOT. DSET(I_LG)%l_vdaer) THEN
          CALL info_bi('NO TRACER IN LAGRANGIAN SPACE with ndrydep=ON',substr)
          CALL info_bi('SWITCHING OF  L_LG', substr)
          L_LG = .FALSE.
       ENDIF
    ENDIF
#endif
!!#D attila -

    IF (.not. L_GP .AND. .not. L_LG .AND. (N_predepvel == 0)) THEN
       DEALLOCATE(DSET)
       CALL info_bi(' NO TRACER AT ALL AVAILABLE FOR DRY DEPOSITION', substr)
       CALL error_bi(' PLEASE SWITCH OFF DRYDEP                    ', substr)
    ENDIF

    IF (p_parallel_io) THEN
       WRITE(*,*) '  Scanning tracers for aerosol models ...'
    END IF

    DO i=1,2
    !------------------------------------------------------------------------
    !   FOR AEROSOL
    !------------------------------------------------------------------------
    ! (2) SCAN FOR RUNNING AEROSOL-MODELS
       IF (i==I_GP) THEN
          IF (.NOT.L_GP) CYCLE
          ti => ti_gp
          DSET(i)%ntrac = ntrac_gp
          DSET(i)%reprid = GP_2D_HORIZONTAL
       ENDIF
       IF (i==I_LG) THEN
          IF (.NOT. L_LG) CYCLE
          ti => ti_lg
          DSET(i)%ntrac = ntrac_lg
          DSET(i)%reprid = LG_ATTILA
       ENDIF

       ALLOCATE(DSET(i)%traermod(DSET(i)%ntrac))
       DSET(i)%traermod(:) = 0
    
       DO jt=1, DSET(i)%ntrac
          ! CHECK IF AEROSOL MODEL
          IF (TRIM(ti(jt)%tp%meta%cask_s(S_aerosol_model)) == '') CYCLE
          lnew = .TRUE.
          ! CHECK IF AEROSOL MODEL ALREADY IN LIST
          DO naermo = 1, DSET(i)%aeromodnum
             IF ( TRIM(ti(jt)%tp%meta%cask_s(S_aerosol_model)) &
                  == DSET(i)%aermodname(naermo) ) THEN
                lnew = .FALSE.
                EXIT
             END IF
          END DO
          ! NEW ENTRY IN LIST
          IF (lnew) THEN
             DSET(i)%aeromodnum = DSET(i)%aeromodnum + 1
             IF (DSET(i)%aeromodnum > NMAX_AEROMODELS) &
                  CALL error_bi( &
                  'RECOMPILATION WITH INCREASED NMAX_AEROMODELS REQUIRED' &
                  , substr)
             DSET(i)%aermodname(DSET(i)%aeromodnum)   = &
                  TRIM(ti(jt)%tp%meta%cask_s(S_aerosol_model))
             DSET(i)%aermodmethod(DSET(i)%aeromodnum) = &
                  ti(jt)%tp%meta%cask_i(I_aerosol_method)
             naermo = DSET(i)%aeromodnum
          END IF
          DSET(i)%traermod(jt) = naermo
       END DO
       
       IF (p_parallel_io) THEN
          WRITE(*,*) '  ... ',DSET(i)%aeromodnum,' aerosol model(s) found:'
          WRITE(*,*) DSET(i)%aermodname(1:DSET(i)%aeromodnum)
       END IF

       ! (3) INITIALIZE MEMORY
       ! (3a) ALLOCATE POINTERS
       ALLOCATE(DSET(i)%radii(DSET(i)%aeromodnum))
       ALLOCATE(DSET(i)%densaer(DSET(i)%aeromodnum))
       ALLOCATE(DSET(i)%sigma(DSET(i)%aeromodnum))
       ALLOCATE(DSET(i)%nmod(DSET(i)%aeromodnum))
       ALLOCATE(DSET(i)%l_anymode(DSET(i)%aeromodnum))
       ALLOCATE(DSET(i)%l_anymass(DSET(i)%aeromodnum))
       ALLOCATE(DSET(i)%l_anynum(DSET(i)%aeromodnum))
       ALLOCATE(DSET(i)%vdrydep_mass(DSET(i)%aeromodnum))
       ALLOCATE(DSET(i)%vdrydep_num(DSET(i)%aeromodnum))
       IF (i==I_LG) THEN
          ALLOCATE(DSET(i)%radii_lg(DSET(i)%aeromodnum))
          ALLOCATE(DSET(i)%densaer_lg(DSET(i)%aeromodnum))
       ENDIF
       ! (3d) SUBMODEL CHANNEL
       IF (DSET(i)%l_vdbl .OR. DSET(i)%l_vdaer) THEN 
          CALL new_channel(status, modstr//'_'//setname(i) &
                          , reprid=DSET(i)%reprid)
          CALL channel_halt(substr, status)
       END IF
       
       !----------------------------------------------------------------------
       !   FOR TRACERS 
       !----------------------------------------------------------------------
       
       if_vdbl: IF (DSET(i)%l_vdbl) THEN
          !----------------------------------------------------------------
          !   checking for which tracer  vdrydep must be calculated 
          !----------------------------------------------------------------
          ! ONLINE BIG-LEAF DRY DEPOSITION VELOCITIES
          
          ! IF TRACERS EXIST, WE NEED TO KNOW WHICH DRY DEPOSITION INDEX
          ! BELONGS TO WHICH TRACER; INITIALIZE THIS FIELD
          IF (DSET(i)%ntrac > 0) THEN
             ALLOCATE (DSET(i)%idt_vdbl(DSET(i)%ntrac))
             DSET(i)%idt_vdbl(:)=0
             ALLOCATE (DSET(i)%idt_diag(DSET(i)%ntrac))
             DSET(i)%idt_diag(:)=0
          END IF
          
          ! INITIALIZE NUMBER OF BIG-LEAF DRY DEPOSITION FIELDS
          nvdbl=0
          ! INITIALIZE NUMBER OF TRACERS FOR DRY DEPOSITION 
          na=0
          ! COUNT DRY DEPOSITION FIELDS AND SAVE INDICES FOR 
          ! CORRESPONDING TRACERS
          ! THE CORRESPONDING TRACER MUST EXIST AND SOME PARAMETERS MUST BE SET
          DO jt = 1, DSET(i)%ntrac
             IF ((ti(jt)%tp%meta%cask_i(I_drydep) == ON)) THEN
                na = na + 1
                DSET(i)%idt_diag(jt) = na
                IF (ti(jt)%tp%ident%medium == AIR) THEN
                   IF ( (ti(jt)%tp%meta%cask_r(R_pss  ) > 0.0) .AND. &
                        (ti(jt)%tp%meta%cask_r(R_molarmass) > 0.0) ) THEN
                      nvdbl = nvdbl + 1
                      ! ADD nvdbl HERE, BECAUSE OF RUNNING INDEX FOR 
                      ! vdbl-fields
                      DSET(i)%idt_vdbl(jt)=nvdbl
                   ELSE IF &
                        (TRIM(ti(jt)%tp%ident%fullname) == 'H2SO4') THEN
                      DSET(i)%idt_h2so4 = jt
                      nvdbl = nvdbl + 1
                      DSET(i)%idt_vdbl(jt)=nvdbl

                      ! a temporary measure !
                   ELSE IF &
                      ( INDEX(TRIM(ti(jt)%tp%ident%fullname),'H2SO4',.TRUE.) == &
                      LEN_TRIM(ti(jt)%tp%ident%fullname)-LEN_TRIM('H2SO4')+1 ) THEN
                      nvdbl = nvdbl + 1
                      DSET(i)%idt_vdbl(jt)=nvdbl

                   ELSE
                   CALL info_bi( &
                        'YOU HAVE NOT DEFINED MOLAR MASS AND/OR HENRY COEFF.'&
                        &//'FOR '//TRIM(ti(jt)%tp%ident%fullname), '  ')
                   ENDIF
                END IF
             ENDIF

#ifdef MESSYTENDENCY
             IF (i==I_GP) THEN
                CALL mtend_register(my_handle, jt)
             ENDIF
#endif

          END DO
          
          DSET(i)%nodoc = na 
          
          ! redefine nodoc and idt_diag if smaller output subset
          ! and l_diagflux=.FALSE.
          
          IF (.not. l_diagflux) THEN
             ! put only named dry deposition fluxes into channels
             DSET(i)%nodoc = 0
             call strcrack(outflux, ';', strname, DSET(i)%nodoc)  
             CALL int2str(teststr, DSET(i)%nodoc, '0', 'X')
             CALL info_bi ('nodoc_'//setname(i)//' : '//teststr , substr)
             DSET(i)%idt_diag(:)=0
             IF (DSET(i)%nodoc /= 0) THEN  
                ! THE CORRESPONDING TRACER MUST EXIST 
                offset = 0 ! setting corrector to zero
                channelloop: DO ic = 1, DSET(i)%nodoc 
                   trac_loop: DO jt = 1, DSET(i)%ntrac
                      IF ((ti(jt)%tp%meta%cask_i(I_drydep) == ON) .AND. &
                           TRIM(ti(jt)%tp%ident%fullname) == &
                           TRIM(strname(ic)))&
                           THEN
                         DSET(i)%idt_diag(jt) = ic - offset
                         CALL int2str(teststr, jt, '0', 'X')
                         CALL info_bi('Tracer '//teststr//&
                              &' identified as'//strname(ic), substr)
                         exit
                      ENDIF
                   ENDDO trac_loop
                   IF (jt > DSET(i)%ntrac) THEN
                      CALL info_bi(' NO corresponding tracer found'&
                           &//'for'//TRIM(strname(ic))//&
                           &' correcting number of channels',substr)
                      offset = offset + 1
                   ENDIF
                ENDDO channelloop
                DSET(i)%nodoc = DSET(i)%nodoc - offset
                CALL int2str(teststr, DSET(i)%nodoc, '0', 'X')
                CALL info_bi(&
                     'nodoc_'//setname(i)//' corrected to: '//teststr, substr)
                IF (DSET(i)%nodoc > 0) l_fluxout = .TRUE.
             ENDIF ! nodoc > 0
          ENDIF ! l_diagflux

          IF (nvdbl > 0) THEN
             ALLOCATE (DSET(i)%vd_bigl(nvdbl))
          ENDIF
          
          !-----------------------------------------------------------------
          !   ASSIGN OUTPUT FIELDS
          !-----------------------------------------------------------------
          
          ! ASSIGNING THE BIG-LEAF DRY DEPOSITION FIELDS
          ! ALLOCATE big-leaf dry deposition channel object pointer
          
          ii = 0  ! big-leaf dry deposition channel objects
          DO jt = 1, DSET(i)%ntrac
             IF (DSET(i)%idt_vdbl(jt) /= 0) THEN
                ii=ii+1
                CALL new_channel_object(status, modstr//'_'//setname(i)  &
                     , 'Vd_'//TRIM(ti(jt)%tp%ident%fullname)       &
                     , p2 = DSET(i)%vd_bigl(ii)%ptr                &
                     , reprid = GP_2D_HORIZONTAL)
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//'_'//setname(i) &
                     , 'Vd_'//TRIM(ti(jt)%tp%ident%fullname)       &
                     , 'units', c='cm s-1')
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//'_'//setname(i) &
                     , 'Vd_'//TRIM(ti(jt)%tp%ident%fullname)       &
                     , 'long_name', c='Big-leaf dry dep. veloc. '&
                     &//TRIM(ti(jt)%tp%ident%fullname))
                CALL channel_halt(substr, status)
             END IF
          END DO
       ENDIF if_vdbl

       ! put all dry deposition fluxes into channels
       IF (l_diagflux .or. l_fluxout) THEN
          ALLOCATE (DSET(i)%drydepflux(DSET(i)%nodoc))
          IF (i == I_LG) ALLOCATE(DSET(i)%drydepflux_lggp(DSET(i)%nodoc))

          DO jt =1, DSET(i)%ntrac
             IF (DSET(i)%idt_diag(jt) /= 0) THEN
                CALL new_channel_object(status, modstr//'_'//setname(i)  &
                     , 'ddepflux_'//TRIM(ti(jt)%tp%ident%fullname)       &
                     , p2 = DSET(i)%drydepflux(DSET(i)%idt_diag(jt))%ptr &
                     , reprid = DSET(i)%reprid & 
                     )
                CALL channel_halt(substr, status)
                SELECT CASE(ti(jt)%tp%ident%quantity)
                CASE(NUMBERDENSITY)
                   unit = '# m-2 s-1'
                CASE(AMOUNTFRACTION)
                   unit = 'molec m-2 s-1'
                CASE DEFAULT
                   unit = 'unknown'
                END SELECT
                CALL new_attribute(status, modstr//'_'//setname(i)  &
                     , 'ddepflux_'//TRIM(ti(jt)%tp%ident%fullname)  &
                     , 'units', c=TRIM(unit)                        &
                     )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//'_'//setname(i)  &
                     , 'ddepflux_'//TRIM(ti(jt)%tp%ident%fullname)  &
                     , 'long_name', c='dry deposition flux for '&
                     &//TRIM(ti(jt)%tp%ident%fullname) )
                CALL channel_halt(substr, status)

                IF (i==I_LG) THEN
                   CALL new_channel_object(status, modstr//'_'//setname(i)    &
                        , 'ddepflux_lggp_'//TRIM(ti(jt)%tp%ident%fullname)    &
                        ,p2=DSET(i)%drydepflux_lggp(DSET(i)%idt_diag(jt))%ptr &
                        , reprid = GP_2D_HORIZONTAL                           &
                        )
                   CALL channel_halt(substr, status)
                   SELECT CASE(ti(jt)%tp%ident%quantity)
                   CASE(NUMBERDENSITY)
                      unit = '# m-2 s-1'
                   CASE(AMOUNTFRACTION)
                      unit = 'molec m-2 s-1'
                   CASE DEFAULT
                      unit = 'unknown'
                   END SELECT
                   CALL new_attribute(status, modstr//'_'//setname(i)      &
                        , 'ddepflux_lggp_'//TRIM(ti(jt)%tp%ident%fullname) &
                        , 'units', c=TRIM(unit)                            &
                        )
                   CALL channel_halt(substr, status)
                   CALL new_attribute(status, modstr//'_'//setname(i)       &
                        , 'ddepflux_lggp_'//TRIM(ti(jt)%tp%ident%fullname)  &
                        , 'long_name', c='dry deposition flux for '&
                        &//TRIM(ti(jt)%tp%ident%fullname) )
                   CALL channel_halt(substr, status)
                   
                ENDIF
             ENDIF
          ENDDO
          
          IF (l_diagflux) THEN 
             ALLOCATE (DSET(i)%drydepfluxsum(DSET(i)%nodoc))
             IF (i==I_LG) ALLOCATE (DSET(i)%drydepfluxsum_lggp(DSET(i)%nodoc))
             DO jt =1, DSET(i)%ntrac
                IF(DSET(i)%idt_diag(jt) /= 0) THEN
                   CALL new_channel_object(status, modstr//'_'//setname(i) &
                        , 'ddepfluxsum_'//TRIM(ti(jt)%tp%ident%fullname)  &
                        ,p2 = DSET(i)%drydepfluxsum(DSET(i)%idt_diag(jt))%ptr &
                        , reprid = DSET(i)%reprid &
                        , lrestreq = .TRUE. )
                   CALL channel_halt(substr, status)
                   SELECT CASE(ti(jt)%tp%ident%quantity)
                   CASE(NUMBERDENSITY)
                      unit = '# m-2'
                   CASE(AMOUNTFRACTION)
                      unit = 'kg(tracer) m-2'
                   CASE DEFAULT
                      unit = 'unknown'
                   END SELECT
                   CALL new_attribute(status, modstr//'_'//setname(i)     &
                        , 'ddepfluxsum_'//TRIM(ti(jt)%tp%ident%fullname)  &
                        , 'units', c=TRIM(unit) )
                   CALL channel_halt(substr, status)
                   CALL new_attribute(status, modstr//'_'//setname(i)     &
                        , 'ddepfluxsum_'//TRIM(ti(jt)%tp%ident%fullname)  &
                        , 'long_name', c='dry deposition flux sum for '&
                        &//TRIM(ti(jt)%tp%ident%fullname) )
                   CALL channel_halt(substr, status)
                   
                   IF (i==I_LG) THEN
                   CALL new_channel_object(status, modstr//'_'//setname(i)    &
                   , 'ddepfluxsum_lggp_'//TRIM(ti(jt)%tp%ident%fullname)      &
                   ,p2 = DSET(i)%drydepfluxsum_lggp(DSET(i)%idt_diag(jt))%ptr &
                   , reprid = GP_2D_HORIZONTAL, lrestreq = .TRUE. )
                   CALL channel_halt(substr, status)
                   SELECT CASE(ti(jt)%tp%ident%quantity)
                   CASE(NUMBERDENSITY)
                      unit = '# m-2'
                   CASE(AMOUNTFRACTION)
                      unit = 'kg(tracer) m-2'
                   CASE DEFAULT
                      unit = 'unknown'
                   END SELECT
                   CALL new_attribute(status, modstr//'_'//setname(i)     &
                   , 'ddepfluxsum_lggp_'//TRIM(ti(jt)%tp%ident%fullname)  &
                   , 'units', c=TRIM(unit) )
                   CALL channel_halt(substr, status)
                   CALL new_attribute(status, modstr//'_'//setname(i)     &
                   , 'ddepfluxsum_lggp_'//TRIM(ti(jt)%tp%ident%fullname)  &
                   , 'long_name', c='dry deposition flux sum for '&
                   &//TRIM(ti(jt)%tp%ident%fullname) )
                   CALL channel_halt(substr, status)
                   ENDIF
                ENDIF
             ENDDO
          END IF
       ENDIF

       ! ALLOCATE MEMORY FOR CORE VARIABLES
       ALLOCATE (DSET(i)%zvdrydep(nproma,ngpblks,DSET(i)%ntrac))
       ALLOCATE (DSET(i)%vd_eff(nproma,ngpblks,DSET(i)%ntrac))
       IF (i==I_LG) THEN
          ALLOCATE(DSET(i)%zvdrydep_lg(NCELL,DSET(i)%ntrac))
          ALLOCATE(DSET(i)%vd_eff_lg(NCELL,DSET(i)%ntrac))
          DSET(i)%zvdrydep_lg(:,:) = 0.
          DSET(i)%vd_eff_lg(:,:)   = 0.
       ENDIF
       DSET(i)%zvdrydep(:,:,:) = 0.
       DSET(i)%vd_eff(:,:,:)   = 0.
       
       ! =====================================================================
       !--- attribute specific parameters in the following order:
       !    - diffusivity coefficient, to correct stomatal resistance for
       !      differences in diffusivity between water vapour and the
       !      specific trace gas (sqrt(molmass trace gas)/sqrt(molmass h2o))
       !    - mesophyll resistance
       !    - cuticle resistance
       !    - soil resistance
       !    - wet skin reservoir resistance
       !    - sea water resistance, which is generally similar to the wet skin
       !    - snow resistance
       ALLOCATE(DSET(i)%diff  (DSET(i)%ntrac)); DSET(i)%diff(:)   = 0.0_dp
       ALLOCATE(DSET(i)%diffrb(DSET(i)%ntrac)); DSET(i)%diffrb(:) = 0.0_dp
       ALLOCATE(DSET(i)%rmes  (DSET(i)%ntrac)); DSET(i)%rmes(:)   = 0.0_dp
       ALLOCATE(DSET(i)%rsoil (DSET(i)%ntrac)); DSET(i)%rsoil(:)  = 0.0_dp
       ALLOCATE(DSET(i)%rwater(DSET(i)%ntrac)); DSET(i)%rwater(:) = 0.0_dp
       ALLOCATE(DSET(i)%rsnow (DSET(i)%ntrac)); DSET(i)%rsnow(:)  = 0.0_dp

       ! ju_te_20180625+ 
       ! variables for output of drydep_calc_rcut
       IF (l_ganzeori) THEN
          ALLOCATE(DSET(i)%rcut  (DSET(i)%ntrac)) ; DSET(i)%rcut(:) = 0.0_dp
          ALLOCATE(DSET(i)%rws   (DSET(i)%ntrac)) ; DSET(i)%rws(:)  = 0.0_dp
       ELSE                                               
          ALLOCATE(DSET(i)%rws_2d  (DSET(i)%ntrac))
          ALLOCATE(DSET(i)%rcut_2d (DSET(i)%ntrac))
          DO jt=1,DSET(i)%ntrac
             IF (DSET(i)%idt_vdbl(jt)/=0) THEN
                CALL new_channel_object(status, modstr//'_'//setname(i) &
                     , 'rws_2d'//TRIM(ti(jt)%tp%ident%fullname)         &
                     , p2 = DSET(i)%rws_2d(jt)%ptr                      &
                     , reprid = GP_2D_HORIZONTAL)
                CALL channel_halt(substr, status)
                CALL new_channel_object(status, modstr//'_'//setname(i) &
                     , 'rcut_2d'//TRIM(ti(jt)%tp%ident%fullname)        &
                     , p2 = DSET(i)%rcut_2d(jt)%ptr                     &
                     , reprid = GP_2D_HORIZONTAL)
                CALL channel_halt(substr, status)
             END IF
          END DO
       END IF
       ! ju_te_20180625-

    ENDDO ! Loop over drydep sets

    ! Allocate memory for surface roughness
    ALLOCATE(z0m(nproma,ngpblks))

    ! predefined dry deposition: diagnostics
    IF (N_predepvel > 0) THEN

       CALL get_channel_info(status, modstr//'_'//setname(1))
       IF (status == 3003) THEN  ! channel does not exist
          CALL new_channel(status, modstr//'_'//setname(1))
       END IF
       CALL channel_halt(substr, status)

       DO jt=1, N_predepvel !loop over number of surface sinks
          IF (X_predepvel(jt)%import_predepvel%TRLOGICAL) THEN

             CALL new_channel_object(status, modstr//'_'//setname(1)  &
                  , 'predef_drydflux_'//&
                  &TRIM(X_predepvel(jt)%import_predepvel%OBJECT)//'_' &
                  &//TRIM(X_predepvel(jt)%import_predepvel%TRNAME)    &
                  , p2 = X_predepvel(jt)%diagflux     &
                  , reprid = GP_2D_HORIZONTAL         &
                  )
             CALL channel_halt(substr, status)
             
             CALL new_attribute(status, modstr//'_'//setname(1) &
                  , 'predef_drydflux_'//&
                  &TRIM(X_predepvel(jt)%import_predepvel%OBJECT)//'_' &
                  &//TRIM(X_predepvel(jt)%import_predepvel%TRNAME)    &
                  , 'units', c='molec m-2 s-1')
             CALL channel_halt(substr, status)
             
             CALL new_channel_object(status, modstr//'_'//setname(1)  &
                  , 'predef_drydsum_'//&
                  &TRIM(X_predepvel(jt)%import_predepvel%OBJECT)//'_' &
                  &//TRIM(X_predepvel(jt)%import_predepvel%TRNAME)    &
                  , p2 = X_predepvel(jt)%diagsum      &
                  , reprid = GP_2D_HORIZONTAL, lrestreq = .TRUE.)
             CALL channel_halt(substr, status)
             
             CALL new_attribute(status, modstr//'_'//setname(1) &
                  , 'predef_drydsum_'//&
                  &TRIM(X_predepvel(jt)%import_predepvel%OBJECT)//'_' &
                  &//TRIM(X_predepvel(jt)%import_predepvel%TRNAME)    &
                  , 'units', c='kg m-2')
             CALL channel_halt(substr, status)
             
          ENDIF

#ifdef MESSYTENDENCY
          CALL get_tracer(status, setname(1) &
               , TRIM(X_predepvel(jt)%name), TRIM(X_predepvel(jt)%subname) &
               , idx=idt)
          CALL tracer_halt(substr,status)
          CALL mtend_register(my_handle, idt)
#endif

       ENDDO
    ENDIF

    CALL end_message_bi(modstr, 'MEMORY INITIALIZATION', substr)

  END SUBROUTINE ddep_init_memory

  !==========================================================================

  SUBROUTINE ddep_init_coupling

    ! BML/MESSy
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, ngpblks, nlev
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_1LEV, LG_ATTILA
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt

    ! MESSy
    USE messy_main_channel,          ONLY: get_channel_info, get_channel_object&
                                         , new_channel_object, new_attribute   &
                                         , STRLEN_OBJECT
    USE messy_main_tracer,           ONLY: AMOUNTFRACTION, NUMBERDENSITY, BIN &
                                         , AEROSOL, t_trinfo_tp               &
                                         , I_TAG_REG_IDT
    USE messy_main_tools,            ONLY: int2str


    IMPLICIT NONE

    INTRINSIC TRIM, SIZE

    ! LOCAL
    ! auxiliary pointer for data managment
    REAL(dp), DIMENSION(:,:,:,:), POINTER ::  mem
    CHARACTER(len=*), PARAMETER  :: substr='ddep_init_coupling'
    INTEGER                      :: status
    INTEGER                      :: jm, jt, zmode, naermo, i, modenum
    INTEGER                      :: jtr ! mz_sg_20100320
    CHARACTER(LEN=STRLEN_OBJECT) :: name
    CHARACTER(LEN=2)             :: modestr
    CHARACTER(LEN=24)            :: trname
    REAL(dp)                     :: henry
    REAL(dp)                     :: molweight
    REAL(dp)                     :: dryreac
    CHARACTER(LEN=STRLEN_MEDIUM+3) :: aermodname_rpr

    TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti => NULL()

    !===================== vdaer ============================

    DO i=1,2
       IF (i==I_GP) THEN
          IF (.NOT.L_GP) CYCLE
          ti => ti_gp
          DSET(i)%reprid = GP_3D_1LEV
       ENDIF
       IF (i==I_LG) THEN
          IF (.NOT. L_LG) CYCLE
          ti => ti_lg
          DSET(i)%reprid =LG_ATTILA
       ENDIF

       CALL start_message_bi(modstr, 'COUPLING '//setname(i), substr)

       if_vdaer: IF (DSET(i)%l_vdaer) THEN
       
          aermo_loop: DO naermo = 1, DSET(i)%aeromodnum

             aermodname_rpr = TRIM(DSET(i)%aermodname(naermo))//'_'//setname(i)
             
             !DSET(i)%l_running(naermo)   = .FALSE.
             DSET(i)%nmod(naermo)        = 0
             DSET(i)%radii(naermo)%ptr   => NULL()
             DSET(i)%densaer(naermo)%ptr => NULL()
             DSET(i)%sigma(naermo)%ptr   => NULL()
             DSET(i)%vdrydep_mass(naermo)%ptr => NULL()
             DSET(i)%vdrydep_num(naermo)%ptr  => NULL()
             IF (i==I_LG) THEN
                DSET(i)%radii_lg(naermo)%ptr   => NULL()
                DSET(i)%densaer_lg(naermo)%ptr => NULL()
             ENDIF
             
             CALL get_channel_info(status, TRIM(aermodname_rpr))
             IF (status /= 0) CALL error_bi(&
              'requested aerosol model '//TRIM(aermodname_rpr)//' not running!'&
               , substr)

            CALL info_bi(&
                  ' calculate dry deposition for '//TRIM(aermodname_rpr) &
                  , substr )
             !
             IF (i==I_GP) THEN
                CALL get_channel_object(status, TRIM(aermodname_rpr) &
                     , 'wetradius', p4=DSET(i)%radii(naermo)%ptr) 
                IF (status /= 0) CALL error_bi( &
                     ' ambient radius not available for '&
                     &//TRIM(aermodname_rpr) , substr)
                !
                CALL get_channel_object(status, TRIM(aermodname_rpr) &
                     , 'densaer', p4=DSET(i)%densaer(naermo)%ptr)
                IF (status /=0 ) CALL error_bi( &
                     ' aerosol density not available for '&
                     &//TRIM(aermodname_rpr) , substr) 
             ELSE IF (i==I_LG) THEN
                CALL get_channel_object(status, TRIM(aermodname_rpr) &
                     , 'wetradius', p2=DSET(i)%radii_lg(naermo)%ptr) 
                IF (status /= 0) CALL error_bi( &
                     ' ambient radius not available for '&
                     &//TRIM(aermodname_rpr)  , substr)
                !
                CALL get_channel_object(status, TRIM(aermodname_rpr) &
                     , 'densaer', p2=DSET(i)%densaer_lg(naermo)%ptr)
                IF (status /=0 ) CALL error_bi( &
                     ' aerosol density not available for '&
                     &//TRIM(aermodname_rpr) , substr) 
                modenum = SIZE(DSET(i)%densaer_lg(naermo)%ptr,2)
                ALLOCATE(DSET(i)%radii(naermo)%ptr(nproma,nlev,modenum,ngpblks))
                ALLOCATE(DSET(i)%densaer(naermo)%ptr(nproma,nlev,modenum,ngpblks))
             ENDIF
             !
          
          CALL get_channel_object(status, TRIM(aermodname_rpr) &
               , 'sigma', p1=DSET(i)%sigma(naermo)%ptr)
          IF (status /=0 ) CALL error_bi( &
               ' sigma not available for '//TRIM(aermodname_rpr) , substr) 
          !
          IF (i==I_GP) THEN
             DSET(i)%nmod(naermo) = SIZE(DSET(i)%radii(naermo)%ptr,_IN_XYZN_)
          ELSE IF (i==I_LG) THEN
             DSET(i)%nmod(naermo) = SIZE(DSET(i)%radii_lg(naermo)%ptr,2)
          ENDIF
          !
          ALLOCATE(DSET(i)%l_anymode(naermo)%lptr(DSET(i)%nmod(naermo)))
          ALLOCATE(DSET(i)%l_anymass(naermo)%lptr(DSET(i)%nmod(naermo)))
          ALLOCATE(DSET(i)%l_anynum(naermo)%lptr(DSET(i)%nmod(naermo)))
          DSET(i)%l_anymode(naermo)%lptr(:) = .FALSE.
          DSET(i)%l_anymass(naermo)%lptr(:) = .FALSE.
          DSET(i)%l_anynum(naermo)%lptr(:)  = .FALSE.
          DO jt = 1, DSET(i)%ntrac
             IF ( (TRIM(ti(jt)%tp%meta%cask_s(S_aerosol_model)) == &
                  TRIM(DSET(i)%aermodname(naermo))) .AND. &
                  (ti(jt)%tp%ident%medium == AEROSOL) ) THEN
                zmode = ti(jt)%tp%meta%cask_i(I_aerosol_mode)
                DSET(i)%l_anymode(naermo)%lptr( zmode ) = .TRUE.
                IF (ti(jt)%tp%ident%quantity == AMOUNTFRACTION) &
                     DSET(i)%l_anymass(naermo)%lptr( zmode ) = .TRUE.
                IF (ti(jt)%tp%ident%quantity == NUMBERDENSITY) &
                     DSET(i)%l_anynum(naermo)%lptr( zmode ) = .TRUE.
             END IF
          END DO
          
          IF (DSET(i)%aermodmethod(naermo) /= BIN) &
               ALLOCATE(DSET(i)%vdrydep_mass(naermo)%ptr(&
               _RI_XYZN_( nproma,ngpblks,1,DSET(i)%nmod(naermo) ) ,1))
          ALLOCATE(DSET(i)%vdrydep_num(naermo)%ptr(&
               _RI_XYZN_(nproma,ngpblks,1,DSET(i)%nmod(naermo)),1))

         mod_loop: DO jm = 1, DSET(i)%nmod(naermo)
             
            CALL int2str(modestr, jm, '0', 'X')

             IF (DSET(i)%aermodmethod(naermo) /= BIN) THEN
                name = &
                  TRIM(DSET(i)%aermodname(naermo))//'_v_mass_m'//TRIM(modestr)
                mem => DSET(i)%vdrydep_mass(naermo)%ptr(_RI_XYZN_(:,:,:,jm),:)
                IF  (DSET(i)%l_anymass(naermo)%lptr(jm)) THEN
                   CALL new_channel_object(status, modstr//'_'//setname(i) &
                        , TRIM(name)                      &
                        , reprid = GP_3D_1LEV             &
                        , mem = mem                       &
                        )
                   CALL channel_halt(substr, status)
                   CALL new_attribute(status, modstr//'_'//setname(i) &
                        , TRIM(name)                 &
                        , 'units', c='m/s' )
                   CALL channel_halt(substr, status)
                END IF
             END IF

             name = TRIM(DSET(i)%aermodname(naermo))//'_v_num_m'//TRIM(modestr)
             mem => DSET(i)%vdrydep_num(naermo)%ptr(_RI_XYZN_(:,:,:,jm),:)
             IF  (DSET(i)%l_anynum(naermo)%lptr(jm)) THEN
                CALL new_channel_object(status, modstr//'_'//setname(i) &
                     , TRIM(name)                      &
                     , reprid = GP_3D_1LEV             &
                     , mem = mem                       &
                     )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//'_'//setname(i) &
                     , TRIM(name)                 &
                     , 'units', c='m/s' )
                CALL channel_halt(substr, status)
             END IF
             
          END DO mod_loop

       END DO aermo_loop
       
    END IF if_vdaer
    ! end aerosol part

    !=START ============== airsea coupling ============================
    CALL info_bi('air-sea coupling in the '//setname(i)//' scheme', modstr)
    CALL get_channel_info(status, 'airsea_'//setname(i))
       IF (status /= 0) THEN
          CALL info_bi(' module AIRSEA not running', substr)  
          CALL info_bi(' no coupling with AIRSEA', substr)
       ELSE
          CALL info_bi(' module AIRSEA running', substr)  
          CALL info_bi(' coupling with AIRSEA', substr)
          ALLOCATE(DSET(i)%l_cpl_airsea(DSET(i)%ntrac))
          DSET(i)%l_cpl_airsea(:) = .FALSE. ! mz_pj_20080804
          DSET(i)%l_airsea_coupling = .FALSE.
          DO jt = 1, DSET(i)%ntrac
             IF (DSET(i)%idt_vdbl(jt) == 0) CYCLE ! mz_pj_20080804
             CALL get_channel_object(status, 'airsea_'//setname(i) &
                  , 'dc_'//TRIM(ti(jt)%tp%ident%basename) )
             IF (status /= 0) THEN
                DSET(i)%l_cpl_airsea(jt) =.FALSE.
             ELSE
                DSET(i)%l_cpl_airsea(jt)  =.TRUE.
                DSET(i)%l_airsea_coupling =.TRUE.
             END IF
          END DO
       END IF
    !=END   ============== airsea coupling ============================


    ! call subroutine to pre-calculate  surface resistances,
    CALL start_message_bi(modstr, 'CALCULATE SURFACE RESISTANCES', substr)

    henry     = 0.
    molweight = 0.
    dryreac   = 0.
    tracer_loop: DO jt = 1, DSET(i)%ntrac
       IF (DSET(i)%idt_vdbl(jt) /= 0) THEN

          ! use name of associated regular tracer for tagged tracers
          jtr = ti(jt)%tp%meta%cask_i(I_TAG_REG_IDT)
          trname    = TRIM(ti(jtr)%tp%ident%basename)
          henry     = ti(jt)%tp%meta%cask_r(R_pss  )
          molweight = ti(jt)%tp%meta%cask_r(R_molarmass)
          dryreac   = ti(jt)%tp%meta%cask_r(R_dryreac_sf)
! ju_te_20180706+
          IF (l_ganzeori) THEN 
             CALL drydep_calc_rs(trname, DSET(i)%ntrac, jt,    &
                  henry, molweight, dryreac, DSET(i)%diff,     &
                  DSET(i)%diffrb, DSET(i)%rmes, DSET(i)%rsoil, &
                  DSET(i)%rwater, DSET(i)%rsnow,               &
                  DSET(i)%rws, DSET(i)%rcut)
          ELSE  
             CALL drydep_calc_rs(trname, DSET(i)%ntrac, jt,    &
                  henry, molweight, dryreac,DSET(i)%diff,      &
                  DSET(i)%diffrb, DSET(i)%rmes, DSET(i)%rsoil, &
                  DSET(i)%rwater, DSET(i)%rsnow )
          END IF
! ju_te_20180706-
       END IF
    END DO tracer_loop
    
    CALL end_message_bi(modstr, 'CALCULATE SURFACE RESISTANCES', substr)

 ENDDO ! set loop

!!#D attila +
#ifdef ECHAM5
    ! Special Lagrangian Information   
    IF (L_LG) THEN
       CALL info_bi('import information from attila submodel ', substr)
       CALL get_channel_object(status, cname='attila', oname='AMCELL' &
                  , p0=MASS_PARCEL )
       CALL channel_halt(substr, status)
       CALL info_bi('import information from attila submodel ', substr)
       CALL get_channel_object(status, cname='attila', oname='PPRESS' &
                  , p1=PRESS_PARCEL )
       CALL channel_halt(substr, status)
    END IF
#endif
!!#D attila -

    !*************************************************************
    ! GET INPUT DATA
    !*************************************************************

    ! get leaf area index
    CALL get_channel_object(status, imp_lai%CHA &
         , imp_lai%OBJ, p2=lai )
    CALL channel_halt(substr, status)

    CALL get_channel_object(status, imp_hc%CHA, &
         imp_hc%OBJ, p2=hc )
    CALL channel_halt(substr, status)

    ! get drag
    CALL get_channel_object(status, imp_drag%CHA &
         , imp_drag%OBJ, p2=drag )
    CALL channel_halt(substr, status)

    ! get soilpH
    CALL get_channel_object(status, imp_soilpH%CHA &
         , imp_soilpH%OBJ, p3=soilpH )
    CALL channel_halt(substr, status)

    IF (.NOT. l_ganzeori) THEN
       CALL get_channel_object(status, rainrate_ls%CHA &
            , rainrate_ls%OBJ, p2=rrain_ls)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, rainrate_cv%CHA &
            , rainrate_cv%OBJ, p2=rrain_cv)
       CALL channel_halt(substr, status)
    END IF

    ! get predefined fields of dry-deposition velocity and affected tracer
    ! NOTE: only implemented for GP-Tracers
    IF ( N_predepvel > 0 ) THEN
       CALL info_bi('predefined dry-deposition velocity fields ... ', substr )
       DO i=1,N_predepvel  !loop over number of surface sinks

          CALL get_channel_object(status, &
               TRIM(X_predepvel(i)%import_predepvel%channel) &
             , TRIM(X_predepvel(i)%import_predepvel%object)  &
             , p2 = X_predepvel(i)%field_predepvel)
          CALL channel_halt(substr, status)

          CALL info_bi('predefined dry-deposition velocity field: '//  &
               TRIM(X_predepvel(i)%import_predepvel%channel) // &
               TRIM(X_predepvel(i)%import_predepvel%object) &
               , substr )

          CALL get_tracer(status, setname(1), TRIM(X_predepvel(i)%name), &
               subname=TRIM(X_predepvel(i)%subname), &
               idx=X_predepvel(i)%trid)
          CALL tracer_halt(substr, status)

          CALL info_bi('predefined-dry-deposition tracer: '//  &
               TRIM(X_predepvel(i)%name) //' '// &
               TRIM(X_predepvel(i)%subname) &
               , substr )

       ENDDO
    ELSE
       CALL info_bi('no predefined fields of dry deposition velocity ', substr )
    ENDIF

    CALL end_message_bi(modstr, 'COUPLING', substr)

 END SUBROUTINE ddep_init_coupling
  
  !==========================================================================
  
  SUBROUTINE ddep_global_start

    ! Author: Astrid Kerkweg, UNI-MZ, 2010

    IMPLICIT NONE

    INTRINSIC EXP, MAX, SQRT

    ! LOCAL
    CHARACTER(len=*), PARAMETER :: substr='ddep_global_start'
    
    z0m(:,:)=100./EXP(1./SQRT(MAX(0.001_dp,drag(:,:))))

  END SUBROUTINE ddep_global_start

  !==========================================================================
  
  SUBROUTINE ddep_vdiff

    USE messy_main_grid_def_mem_bi, ONLY: nlev, nproma, jrow, kproma
    USE messy_main_grid_def_bi,     ONLY: deltaz
    USE messy_main_data_bi,         ONLY: pressi_3d
#if defined (ECHAM5)
    USE messy_main_data_bi,         ONLY: pxtems
#endif
     USE messy_main_timer,        ONLY: time_step_len, delta_time
#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, pxtm1 => qxtm1 & 
                                      , ntrac_gp, ti_gp
#else
    USE messy_main_tracer_mem_bi, ONLY: ntrac_gp, ti_gp
#endif
!!#D attila +
#ifdef ECHAM5
    USE messy_attila_tools_e5,    ONLY: lg2gp_e5, LG2GP_AVE
#endif
!!#D attila -

    ! MESSy
    USE messy_main_tracer,        ONLY: NUMBERDENSITY, AMOUNTFRACTION &
                                      , I_TAG_REG_IDT
    
    CHARACTER(LEN=*), PARAMETER :: substr='ddep_vdiff'

    REAL(dp) :: zdz(nproma)
    REAL(dp) :: zdp(nproma)
#ifdef MESSYTENDENCY
    REAL(dp) :: xtstart(_RI_X_ZN_(nproma,nlev,ntrac_gp))
    REAL(dp) :: zxtte(kproma,nlev)
#endif
    REAL(dp) :: xtp1(kproma,ntrac_gp)
    REAL(dp) :: conv
    REAL(dp) :: zdrydepflux(nproma,DSET(I_GP)%ntrac)
    INTEGER  :: jt, naer
    INTEGER  :: jtr

#if defined (ECHAM5) || defined(CESM1)
    REAL(dp), POINTER :: zxtems(:,:)
#endif

    REAL(dp) :: conv_sum

    gridpoint: IF (L_GP) THEN

    IF (DSET(I_GP)%l_vdbl) THEN
       CALL calc_vd_drydep(I_GP, kproma, jrow)
    ENDIF

#if defined(ECHAM5)
    zxtems => pxtems(:,1,:,jrow)
#endif

#ifndef MESSYTENDENCY
    xtp1(1:kproma,:) = pxtm1(_RI_X_ZN_(1:kproma,nlev,:)) + &
         pxtte(_RI_X_ZN_(1:kproma,nlev,:)) *time_step_len
#else
    CALL mtend_get_start_l(mtend_id_tracer, v0t=xtstart)
    xtp1(1:kproma,:) = xtstart(_RI_X_ZN_(1:kproma,nlev,:))
#endif

    zdrydepflux(:,:) = 0._dp
    zdz(:)           = 0._dp       
    zdp(:)           = 0._dp

    ! CALCULATION OF DRY DEPOSITION FLUXES
    zdz(1:kproma) = deltaz(_RI_XYZ__(1:kproma,jrow,nlev))
    zdp(1:kproma) = pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev+1)) - &
         pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev))

    ! The tracer flux to the surface is  equivalent to the deposition velocity:
    ! The dry deposition flux is calculated from the surface layer concentration
    ! and the dry deposition velocity, which includes the impact of the 
    ! sub-timestep scale decrease in the concentrations due to dry deposition!
    ! This prevents the dry deposition flux to remove more mass then present
    ! in the surface layer.
             
    tracer_loop: DO jt=1,DSET(I_GP)%ntrac
       if_drydep: IF(ti_gp(jt)%tp%meta%cask_i(I_drydep) == ON) THEN
          ! calculate effective deposition velocity
          
          ! to get the units as needed for xtems in vdiff the concentration is
          ! multiplied by vd_eff and density of air 
          ! -> mol/mol  Kg/(m2s) for tracer
          ! ->   1/mol  Kg/(m2s) for particle number
          zdrydepflux(1:kproma,jt)=(xtp1(1:kproma,jt)  &
#if defined(ECHAM5)
               + zxtems(1:kproma,jt)/ (zdp(1:kproma)/g)*time_step_len &
#endif
               ) * zdp(1:kproma)/g /zdz(1:kproma)* &
               DSET(I_GP)%vd_eff(1:kproma,jrow,jt)

          SELECT CASE(ti_gp(jt)%tp%ident%quantity)
          CASE(AMOUNTFRACTION)
             conv = conv_mr
             ! conv_sum = (ti_gp(jt)%tp%meta%cask_r(R_molarmass)/M_air) ! old
             ! use molar mass of associated regular tracer for tagged tracers
             jtr = ti_gp(jt)%tp%meta%cask_i(I_TAG_REG_IDT)
             conv_sum = (ti_gp(jtr)%tp%meta%cask_r(R_molarmass)/M_air)
          CASE(NUMBERDENSITY)
             conv = conv_nu
             conv_sum = 1.0E-03/M_air
          CASE DEFAULT
             conv = 0.0_dp
             conv_sum = 0.0_dp
          END SELECT
          
          IF (l_diagflux) THEN
             DSET(I_GP)%drydepflux(DSET(I_GP)%idt_diag(jt))%ptr(1:kproma,jrow)=&
                  zdrydepflux(1:kproma,jt) * conv
             DSET(I_GP)%drydepfluxsum(&
                    DSET(I_GP)%idt_diag(jt))%ptr(1:kproma,jrow) = &
                  DSET(I_GP)%drydepfluxsum(&
                    DSET(I_GP)%idt_diag(jt))%ptr(1:kproma,jrow)&
                  + zdrydepflux(1:kproma,jt) &
                  * delta_time * conv_sum
          ELSEIF (l_fluxout) THEN
             IF (DSET(I_GP)%idt_diag(jt) /= 0) &
                  DSET(I_GP)%drydepflux(&
                     DSET(I_GP)%idt_diag(jt))%ptr(1:kproma,jrow)= &
                  zdrydepflux(1:kproma,jt) * conv
          END IF

#if defined(ECHAM5)
          IF (.NOT. l_tendency ) THEN
             ! units mol (trac.) mol-1 (air) kg(air) m-2 s-1
             zxtems(1:kproma,jt)=zxtems(1:kproma,jt) &
                  - zdrydepflux(1:kproma,jt)
          ELSE
#endif
#ifndef MESSYTENDENCY
             pxtte(_RI_X_ZN_(1:kproma,nlev,jt)) = pxtte(_RI_X_ZN_(1:kproma,nlev,jt))      &
                  - zdrydepflux(1:kproma,jt)/ (zdp(1:kproma)/g)
#else
             zxtte(:,:) = 0._dp
             zxtte(1:kproma,nlev) = &
                  - zdrydepflux(1:kproma,jt)/ (zdp(1:kproma)/g)
             CALL mtend_add_l(my_handle, jt, px=zxtte)
#endif
#if defined (ECHAM5)
          END IF  ! .not.l_tendency
#endif          

       END IF if_drydep
    END DO tracer_loop

    END IF gridpoint

!!#D attila +
#ifdef ECHAM5
    lagrange: IF (L_LG) THEN

    first_row: IF (jrow == 1) THEN

    DO naer=1,DSET(I_LG)%aeromodnum
       call lg2gp_e5(DSET(I_LG)%radii_lg(naer)%ptr,&
            DSET(I_LG)%radii(naer)%ptr, &
            LG2GP_AVE, fill_value=0._dp)
       call lg2gp_e5(DSET(I_LG)%densaer_lg(naer)%ptr,&
            DSET(I_LG)%densaer(naer)%ptr, &
            LG2GP_AVE, fill_value=0._dp)
    ENDDO

    END IF first_row

    IF (DSET(I_LG)%l_vdbl) THEN
          CALL calc_vd_drydep(I_LG, kproma, jrow)
    ENDIF

    END IF lagrange
#endif
!!#D attila -

    IF (N_predepvel > 0) THEN  !number of surface sinks
       CALL ddep_predef_vdiff
    ENDIF

  END SUBROUTINE ddep_vdiff

  !==========================================================================
  SUBROUTINE ddep_global_end
!!#D attila +
#ifdef ECHAM5
    USE messy_main_grid_def_mem_bi, ONLY: nlev, nproma, npromz, ngpblks
    USE messy_main_grid_def_bi,     ONLY: gboxarea_2d, deltaz
    USE messy_main_data_bi,         ONLY: pressi_3d
    USE messy_main_timer,         ONLY: time_step_len, delta_time       
    USE messy_main_tracer_mem_bi, ONLY: pxtte_a => qxtte_a,pxtm1_a => qxtm1_a & 
                                      , ntrac_lg, ti_lg, NCELL
    USE messy_attila_tools_e5,    ONLY: gp2lg_e5, lg2gp_e5, LG2GP_AVE,        &
                                         lggpte2lgte_e5
    USE messy_main_tracer,        ONLY: AEROSOL, NUMBERDENSITY, AMOUNTFRACTION &
                                      , I_TAG_REG_IDT

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, MAX

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='ddep_global_end'
    REAL(dp) :: conv
    REAL(dp), DIMENSION(:,:), POINTER :: zdz         => NULL()
    REAL(dp), DIMENSION(:,:), POINTER :: zdp         => NULL()
    REAL(DP), DIMENSION(:),   POINTER :: zdz_lg      => NULL()
    REAL(DP), DIMENSION(:),   POINTER :: zdp_lg      => NULL()
    REAL(DP), DIMENSION(:),   POINTER :: zps_lg      => NULL()
    REAL(DP), DIMENSION(:),   POINTER :: gboxarea_lg => NULL() 
    INTEGER  :: jt, jn
    INTEGER  :: jtr
    REAL(dp), DIMENSION(:),     POINTER :: ptr_lg       => NULL()
    REAL(dp), DIMENSION(:),     POINTER :: ptr_lg_te    => NULL()
    REAL(dp), DIMENSION(:,:),   POINTER :: ptr_gp       => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: ptr_3d_gp    => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: ptr_3d_gp_te => NULL()
    ! FOR LAGRANGIAN
    LOGICAL, DIMENSION(:), POINTER :: liplev => NULL()
    REAL(dp), TARGET :: xtp1_lg(NCELL, ntrac_lg)
    REAL(dp), TARGET :: xtp1_lggp(nproma,nlev,ntrac_lg,ngpblks)
    REAL(dp), TARGET :: pxtte_a_lggp(nproma,nlev,ntrac_lg,ngpblks)
    REAL(dp)         :: pxtte_a_lg(NCELL, ntrac_lg)
    REAL(dp)         :: zdrydepflux_lg(NCELL,ntrac_lg)
    REAL(dp)         :: zdrydepflux_lggp(nproma,nlev,ntrac_lg,ngpblks) 
    REAL(dp)         :: ztmst, zdtime
    REAL(dp)         :: conv_sum

    IF (.NOT. L_LG) RETURN

    ztmst  = time_step_len
    zdtime = delta_time
    zdrydepflux_lg(:,:)  =0._dp 
    zdrydepflux_lggp(:,:,:,:)=0._dp 
    ALLOCATE(zdz(nproma,ngpblks))
    ALLOCATE(zdp(nproma,ngpblks))
    ALLOCATE(zdz_lg(NCELL))
    ALLOCATE(zdp_lg(NCELL))
    ALLOCATE(zps_lg(NCELL))
    ALLOCATE(gboxarea_lg(NCELL))
    zdz(:,:) = 0._dp       
    zdp(:,:) = 0._dp
    zdz_lg(:) = 0._dp       
    zdp_lg(:) = 0._dp
    zps_lg(:) = 0._dp       
    gboxarea_lg(:) = 0._dp

    ! CALCULATION OF DRX DEPOSITION FLUXES
    zdz(:,:) = deltaz(:,nlev,:)
    zdp(:,:) = pressi_3d(:,nlev+1,:) - pressi_3d(:,nlev,:)

    xtp1_lg(1:NCELL,:) = pxtm1_a(1:NCELL,:) + pxtte_a(1:NCELL,:) *ztmst

    CALL gp2lg_e5( zdp, zdp_lg)
    CALL gp2lg_e5( zdz, zdz_lg)

    tracer_loop: DO jt = 1, DSET(I_LG)%ntrac

       IF (ti_lg(jt)%tp%meta%cask_i(I_drydep) /= ON) CYCLE
       IF ( (ti_lg(jt)%tp%ident%medium /= AEROSOL) .AND. &
            (DSET(I_LG)%idt_vdbl(jt) == 0) ) CYCLE

       ptr_gp => DSET(I_LG)%zvdrydep(:,:,jt)
       ptr_lg => DSET(I_LG)%zvdrydep_lg(:,jt)
       CALL gp2lg_e5(ptr_gp,ptr_lg, llev=liplev)
                   
       ptr_gp =>  DSET(I_LG)%vd_eff(:,:,jt)
       ptr_lg =>  DSET(I_LG)%vd_eff_lg(:,jt)
       CALL gp2lg_e5( ptr_gp, ptr_lg)

       SELECT CASE(ti_lg(jt)%tp%ident%quantity)
       CASE(AMOUNTFRACTION)
          conv = conv_mr

          ! conv_sum = (ti_lg(jt)%tp%meta%cask_r(R_molarmass)/M_air) ! old
          ! use molar mass of associated regular tracer for tagged tracers
          jtr = ti_lg(jt)%tp%meta%cask_i(I_TAG_REG_IDT)
          conv_sum = (ti_lg(jtr)%tp%meta%cask_r(R_molarmass)/M_air)
       CASE(NUMBERDENSITY)
          conv = conv_nu
          conv_sum = 1.0E-03/M_air
       CASE DEFAULT
          conv = 0.0_dp
          conv_sum = 0.0_dp
       END SELECT

       ! we have four options for the lagrangian deposition
       ! FLUXES
       ! we do it two times to avoid zeroes in the rest of the array...
       ! zxtems(:,jt)/ (zdp(1:kproma)/g)*ztmst)*& !skipped zxtems term!!!!
       ! zxtems(:,jt)/ (zdp(1:kproma)/g)*ztmst)*& !skipped zxtems term!!!!
       SELECT CASE(i_lg_method)
       CASE(1)  
          DO jn=1,NCELL
             IF (liplev(jn)) THEN
                ! mol/mol kg/m2s -->
                ! we multiplied for the density of the grid/parcel 
                ! (we suppouse they are the same!) 
                zdrydepflux_lg(jn,jt) = xtp1_lg(jn,jt)*  &
                     (DSET(I_LG)%vd_eff_lg(jn,jt))/g* &   
                     zdp_lg(jn)/zdz_lg(jn)
                pxtte_a(jn,jt) =  pxtte_a(jn,jt) - &
                     zdrydepflux_lg(jn,jt)/(zdp_lg(jn)/g)
             ELSE 
                zdrydepflux_lg(jn,jt) = 0._dp
                pxtte_a(jn,jt) = pxtte_a(jn,jt)          
             END IF
          END DO
          !-----------------------------------------
       CASE (2)
          CALL gp2lg_e5(gboxarea_2d, gboxarea_lg)
          DO jn=1,NCELL
             IF (liplev(jn)) THEN
                ! mol/mol kg/m2s -->
                ! we multiplied for the density of the grid/parcel 
                ! (we suppouse they are the same!) 
                zdrydepflux_lg(jn,jt) = xtp1_lg(jn,jt)*  &
                     DSET(I_LG)%vd_eff_lg(jn,jt)/g* &   
                     zdp_lg(jn)/zdz_lg(jn)
                pxtte_a(jn,jt) = pxtte_a(jn,jt)          &
                     - zdrydepflux_lg(jn,jt)             &
                     *gboxarea_lg(jn)/MASS_PARCEL 
             ELSE 
                zdrydepflux_lg(jn,jt) = 0._dp
                pxtte_a(jn,jt) = pxtte_a(jn,jt)          
             END IF
          END DO
          !-----------------------------------------
       CASE(3)
          ptr_gp => pressi_3d(:,nlev+1,:)
          CALL gp2lg_e5(ptr_gp, zps_lg)
          !
          ! here we suppose that any parcel has the same density
          ! of the grid box where they are located
          DO jn=1,NCELL
             !here zdp_lg and zdz_lg are still equal to GP!!!!!
             IF (liplev(jn)) THEN
                zdz_lg(jn)=(MAX((-1._dp*(PRESS_PARCEL(jn) - &
                     zps_lg(jn))), 1.E-20_dp) /zdp_lg(jn))*zdz_lg(jn)
                zdp_lg(jn)=(MAX((-1._dp*(PRESS_PARCEL(jn) - &
                     zps_lg(jn))), 1.E-20_dp))      
                DSET(I_LG)%vd_eff_lg(jn,jt) = &
                     drydep_calc_vd_eff(DSET(I_LG)%zvdrydep_lg(&
                     jn,jt),ztmst, zdz_lg(jn))
                ! mol/mol kg/m2s --> 
                ! we multiplied for the density of the grid/parcel 
                ! (we suppouse they are the same!) 
                zdrydepflux_lg(jn,jt) = xtp1_lg(jn,jt)*  &
                     (DSET(I_LG)%vd_eff_lg(jn,jt))/g*  &  
                     zdp_lg(jn)/zdz_lg(jn)
                pxtte_a(jn,jt) = pxtte_a(jn,jt)          &
                     - zdrydepflux_lg(jn,jt) /(zdp_lg(jn)/g)
             ELSE 
                zdrydepflux_lg(jn,jt) = 0._dp
                pxtte_a(jn,jt) = pxtte_a(jn,jt)          
             END IF
          END DO
          !-----------------------------------------
       CASE (4)
          ptr_lg => xtp1_lg(1:NCELL,jt)
          ptr_3d_gp => xtp1_lggp(1:nproma,1:nlev,jt,1:ngpblks)
          CALL lg2gp_e5(ptr_lg, ptr_3d_gp, LG2GP_AVE, fill_value=0._dp)

          zdrydepflux_lggp(1:nproma,nlev,jt,1:ngpblks-1)= &
               xtp1_lggp(1:nproma,nlev,jt,1:ngpblks-1)*  &
               zdp(1:nproma,1:ngpblks-1)/g/zdz(1:nproma,1:ngpblks-1) &
               *DSET(I_LG)%vd_eff(1:nproma,1:ngpblks-1,jt)

          zdrydepflux_lggp(1:npromz,nlev,jt,ngpblks)= &
               xtp1_lggp(1:npromz,nlev,jt,ngpblks)*  &
               zdp(1:npromz,ngpblks)/g/zdz(1:npromz,ngpblks)* &
               DSET(I_LG)%vd_eff(1:npromz,ngpblks,jt)
          ! TRACER TENDENCY
          ! we do it two times to avoid zeroes in the rest of the array...
          pxtte_a_lggp(1:nproma,1:nlev,jt,1:ngpblks) = 0.
          pxtte_a_lggp(1:nproma,nlev,jt,1:ngpblks-1) = &
               - zdrydepflux_lggp(1:nproma,nlev,jt,1:ngpblks-1) &
               / (zdp(1:nproma,1:ngpblks-1)/g)

          pxtte_a_lggp(1:npromz,nlev,jt,ngpblks) = &
               - zdrydepflux_lggp(1:npromz,nlev,jt,ngpblks)   &
               / (zdp(1:npromz,ngpblks)/g)
          !!!!!!!!!
          ptr_lg => xtp1_lg(1:NCELL,jt)
          !ptr_lg_te=> pxtte_a_lg(1:NCELL,jt)
          ptr_3d_gp_te => pxtte_a_lggp(1:nproma,1:nlev,jt,1:ngpblks)
          ptr_3d_gp => xtp1_lggp(1:nproma,1:nlev,jt,1:ngpblks)
          CALL lggpte2lgte_e5(ptr_3d_gp, ptr_3d_gp_te, ptr_lg, ptr_lg_te)
          !                           cell_loop: DO jn=1,NCELL
          !is_level: IF (liplev(jn)) THEN
          pxtte_a_lg(1:NCELL,jt)= ptr_lg_te(1:NCELL)
          pxtte_a(1:NCELL,jt) = pxtte_a(1:NCELL,jt) + pxtte_a_lg(1:NCELL,jt)
          ! mol/mol kg/m2s -->
          IF (l_diagflux) THEN   !output
             DSET(I_LG)%drydepflux_lggp(&
                    DSET(I_LG)%idt_diag(jt))%ptr(1:nproma,1:ngpblks)= &
                  zdrydepflux_lggp(1:nproma,nlev,jt,1:ngpblks)* conv
             DSET(I_LG)%drydepfluxsum_lggp(&
                    DSET(I_LG)%idt_diag(jt))%ptr(1:nproma,1:ngpblks) =  & 
                  DSET(I_LG)%drydepfluxsum_lggp(&
                    DSET(I_LG)%idt_diag(jt))%ptr(1:nproma,1:ngpblks) &
                  + zdrydepflux_lggp(1:nproma,nlev,jt,1:ngpblks) * zdtime     &
                  * conv_sum
          ELSEIF (l_fluxout) THEN
             IF (DSET(I_LG)%idt_diag(jt) /= 0) &
                  DSET(I_LG)%drydepflux_lggp(&
                  DSET(I_LG)%idt_diag(jt))%ptr(:,:) = &
                  zdrydepflux_lggp(:,nlev,jt,:)*conv
          END IF
       CASE DEFAULT
          CALL error_bi('NO VALID METHOD SELECTED FOR LAGRANGIAN'//&
               &' DRY DEPOSITION' , substr)
       END SELECT

       IF (i_lg_method .ne. 4) THEN
          ! final output for analysis is done here.
          ! outside cell loop but inside tracer loop!
          IF (l_diagflux) THEN   !output
             DSET(I_LG)%drydepflux(DSET(I_LG)%idt_diag(jt))%ptr(1:NCELL,1) = &
                  zdrydepflux_lg(1:NCELL,jt) * conv
             ALLOCATE (ptr_3d_gp(1:nproma,1:nlev,1:ngpblks))
             ptr_lg => DSET(I_LG)%drydepflux(&
                    DSET(I_LG)%idt_diag(jt))%ptr(1:NCELL,1)
             CALL lg2gp_e5(ptr_lg,ptr_3d_gp, LG2GP_AVE, fill_value=0._dp)
             DSET(I_LG)%drydepflux_lggp(&
                    DSET(I_LG)%idt_diag(jt))%ptr(1:nproma,1:ngpblks)= &
                  ptr_3d_gp(1:nproma,nlev,1:ngpblks)
             DSET(I_LG)%drydepfluxsum(&
                    DSET(I_LG)%idt_diag(jt))%ptr(1:NCELL,1) =  & 
                  DSET(I_LG)%drydepfluxsum(&
                    DSET(I_LG)%idt_diag(jt))%ptr(1:NCELL,1) &
                    + zdrydepflux_lg(1:NCELL,jt) * zdtime          &
                    * conv_sum
             ptr_lg => DSET(I_LG)%drydepfluxsum(&
                    DSET(I_LG)%idt_diag(jt))%ptr(1:NCELL,1)
             CALL lg2gp_e5(ptr_lg, ptr_3d_gp, LG2GP_AVE, fill_value=0._dp)
             DSET(I_LG)%drydepfluxsum_lggp(&
                  DSET(I_LG)%idt_diag(jt))%ptr(1:nproma,1:ngpblks)= &
                  ptr_3d_gp(1:nproma,nlev,1:ngpblks)
             DEALLOCATE  (ptr_3d_gp)
             NULLIFY(ptr_3d_gp)
          ELSEIF (l_fluxout) THEN
             IF (DSET(I_LG)%idt_diag(jt) /= 0) THEN 
                DSET(I_LG)%drydepflux(&
                  DSET(I_LG)%idt_diag(jt))%ptr(1:NCELL,1) = &
                  zdrydepflux_lg(1:NCELL,jt) * conv
                ptr_lg => DSET(I_LG)%drydepflux(&
                        DSET(I_LG)%idt_diag(jt))%ptr(1:NCELL,1)
                ALLOCATE(ptr_3d_gp(1:nproma,1:nlev,1:ngpblks))
                CALL lg2gp_e5(ptr_lg, ptr_3d_gp, LG2GP_AVE, fill_value=0._dp)
                DSET(I_LG)%drydepflux_lggp(&
                       DSET(I_LG)%idt_diag(jt))%ptr(1:nproma,1:ngpblks)= &
                     ptr_3d_gp(1:nproma,nlev,1:ngpblks)
                DEALLOCATE (ptr_3d_gp)
                NULLIFY(ptr_3d_gp)
             END IF
          END IF
       ENDIF

       ! CLEAN UP
       IF (ASSOCIATED(liplev)) DEALLOCATE(liplev)
       NULLIFY(liplev)

    END DO tracer_loop

    DEALLOCATE(zdz, zdz_lg, zdp, zdp_lg, zps_lg, gboxarea_lg)
#endif
!!#D attila -
  END SUBROUTINE ddep_global_end

  !==========================================================================

  SUBROUTINE ddep_free_memory

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE

    INTEGER :: i, naermo

    DEALLOCATE(z0m) 

    IF (ASSOCIATED(DSET)) THEN
       DO i=1,SIZE(DSET)
          IF (ASSOCIATED(DSET(i)%l_cpl_airsea)) &
               DEALLOCATE(DSET(i)%l_cpl_airsea)
          Do naermo = 1, DSET(i)%aeromodnum
             IF (ASSOCIATED(DSET(i)%vdrydep_num(naermo)%ptr)) &
                  DEALLOCATE(DSET(i)%vdrydep_num(naermo)%ptr)
             IF (ASSOCIATED(DSET(i)%vdrydep_mass(naermo)%ptr)) &
                  DEALLOCATE(DSET(i)%vdrydep_mass(naermo)%ptr)
             IF (ASSOCIATED(DSET(i)%l_anymode(naermo)%lptr)) &
                  DEALLOCATE(DSET(i)%l_anymode(naermo)%lptr)  
             IF (ASSOCIATED(DSET(i)%l_anymass(naermo)%lptr)) &
                  DEALLOCATE(DSET(i)%l_anymass(naermo)%lptr)
             IF (ASSOCIATED(DSET(i)%l_anynum(naermo)%lptr)) &
                  DEALLOCATE(DSET(i)%l_anynum(naermo)%lptr)
          ENDDO
          IF (ASSOCIATED(DSET(i)%diff))     DEALLOCATE(DSET(i)%diff) 
          IF (ASSOCIATED(DSET(i)%diffrb))   DEALLOCATE(DSET(i)%diffrb)
          IF (ASSOCIATED(DSET(i)%rmes))     DEALLOCATE(DSET(i)%rmes)
          IF (ASSOCIATED(DSET(i)%rcut))     DEALLOCATE(DSET(i)%rcut)
          IF (ASSOCIATED(DSET(i)%rsoil))    DEALLOCATE(DSET(i)%rsoil)
          IF (ASSOCIATED(DSET(i)%rws))      DEALLOCATE(DSET(i)%rws)
          IF (ASSOCIATED(DSET(i)%rwater))   DEALLOCATE(DSET(i)%rwater)
          IF (ASSOCIATED(DSET(i)%rsnow))    DEALLOCATE(DSET(i)%rsnow)
          IF (ASSOCIATED(DSET(i)%idt_vdbl)) DEALLOCATE(DSET(i)%idt_vdbl)
          IF (ASSOCIATED(DSET(i)%idt_diag)) DEALLOCATE(DSET(i)%idt_diag)
          IF (ASSOCIATED(DSET(i)%vd_bigl))  DEALLOCATE(DSET(i)%vd_bigl)
          ! ju_te_20180705+
          IF (ASSOCIATED(DSET(i)%rcut_2d))  DEALLOCATE(DSET(i)%rcut_2d)
          IF (ASSOCIATED(DSET(i)%rws_2d))   DEALLOCATE(DSET(i)%rws_2d)
          ! ju_te_20180705-

          IF (ASSOCIATED(DSET(i)%drydepflux)) DEALLOCATE (DSET(i)%drydepflux)
          IF (ASSOCIATED(DSET(i)%drydepfluxsum)) &
               DEALLOCATE (DSET(i)%drydepfluxsum)
          IF (ASSOCIATED(DSET(i)%zvdrydep)) DEALLOCATE (DSET(i)%zvdrydep)
          IF (ASSOCIATED(DSET(i)%vd_eff))   DEALLOCATE (DSET(i)%vd_eff)
          IF (i==I_LG) THEN
             IF (ASSOCIATED(DSET(i)%drydepfluxsum_lggp)) &
                  DEALLOCATE(DSET(i)%drydepfluxsum_lggp)
             IF (ASSOCIATED(DSET(i)%zvdrydep_lg)) &
                  DEALLOCATE(DSET(i)%zvdrydep_lg)
             IF (ASSOCIATED(DSET(i)%vd_eff_lg)) &
                  DEALLOCATE(DSET(i)%vd_eff_lg)
             IF (ASSOCIATED(DSET(i)%radii_lg)) &
                  DEALLOCATE(DSET(i)%radii_lg)
             IF (ASSOCIATED(DSET(i)%densaer_lg)) &
                  DEALLOCATE(DSET(i)%densaer_lg)
          ENDIF
          IF (ASSOCIATED(DSET(i)%radii)) &
               DEALLOCATE(DSET(i)%radii)
          IF (ASSOCIATED(DSET(i)%densaer)) &
               DEALLOCATE(DSET(i)%densaer)
          IF (ASSOCIATED(DSET(i)%sigma))&
               DEALLOCATE(DSET(i)%sigma)
          IF (ASSOCIATED(DSET(i)%nmod)) &
               DEALLOCATE(DSET(i)%nmod)
          IF (ASSOCIATED(DSET(i)%l_anymode)) &
               DEALLOCATE(DSET(i)%l_anymode)
          IF (ASSOCIATED(DSET(i)%l_anymass)) &
               DEALLOCATE(DSET(i)%l_anymass)
          IF (ASSOCIATED(DSET(i)%l_anynum)) &
               DEALLOCATE(DSET(i)%l_anynum)
          IF (ASSOCIATED(DSET(i)%vdrydep_mass)) &
               DEALLOCATE(DSET(i)%vdrydep_mass)
          IF (ASSOCIATED(DSET(i)%vdrydep_num)) &
               DEALLOCATE(DSET(i)%vdrydep_num)
          IF (ASSOCIATED(DSET(i)%traermod)) DEALLOCATE(DSET(i)%traermod)
    ENDDO
       IF (ASSOCIATED(DSET)) DEALLOCATE(DSET)
    ENDIF
  END SUBROUTINE ddep_free_memory


! ===========================================================================
! PRIVATE ROUTINES
! ===========================================================================
  !==========================================================================

  SUBROUTINE calc_vd_drydep(flag, locproma, locrow)

    ! Purpose:
    ! ---------
    ! This subroutine calculates the dry deposition
    ! of aerosols and gases.
    !
    ! Authors:
    ! ----------
    ! Philip Stier,      MPI-MET
    ! Laurens Ganzeveld, MPI-CHEM
    ! Astrid Kerkweg,    MPI-CHEM, UNi-MAINZ
    !
    ! Interface:
    ! ----------
    ! *drydep_calc_drydep* is called from *drydep_vdiff*
    ! or *drydep_global_end"
    !
    ! Method:
    ! -------
    !@@@ Careful - Laurens seems to assume Vd in cm/s!
    !@@@ This is not taken into account!
    
    USE messy_main_grid_def_mem_bi, ONLY: nlev
    USE messy_main_grid_def_bi,     ONLY: deltaz
    USE messy_main_data_bi,         ONLY: cvs, cvw, alake                 &
                                        , seaice, vgrat                   &
                                        , cdnl, cfncl,cfml, az0           &
                                        , um1,  vm1, tvir, tvl, ril       &
                                        , cdnw, cfmw, cfncw, riw, tvw     &
                                        , cdni, cfmi,cfnci, rii, tvi      &
                                        , tslm1, rh_2m                    &
                                        , rco_leaf, fws                   &
                                        , u10, v10                        &
                                        , slf
    USE messy_main_timer,         ONLY: time_step_len 
    USE messy_main_constants_mem, ONLY: g
    USE messy_main_tracer,        ONLY: BIN, t_trinfo_tp                       &
                                      , AEROSOL, NUMBERDENSITY, AMOUNTFRACTION &
                                      , I_TAG_REG_IDT 

    IMPLICIT NONE

    INTRINSIC :: TINY

    INTEGER, INTENT(IN)     :: flag
    INTEGER, INTENT(IN)     :: LOCROW
    INTEGER, INTENT(IN)     :: locproma
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='calc_vd_drydep'
    TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti 
    LOGICAL  :: loland(locproma)
    REAL(dp) :: zrahl(locproma), zrahw(locproma),    zrahi(locproma),   &
         zrahveg(locproma),      zrahslsn(locproma),                    &
         zustarl(locproma),      zustarw(locproma),  zustari(locproma), &
         zustveg(locproma),      zustslsn(locproma),                    &
         zevap(locproma),        zvgrat(locproma),   zcvbs(locproma),   &
         zfrw(locproma),         zfri(locproma)
    REAL(dp) :: zdz(locproma)
    REAL(dp) :: zzvdrydep(locproma)
    REAL(dp) :: abswind(locproma), abswind10(locproma)
    REAL(dp) :: radius_mass(locproma), radius_number(locproma)
    REAL(dp) :: zalpha(locproma) 
    REAL(dp) :: zbubble(locproma)
    REAL(dp) :: zbeta(locproma)
    REAL(dp) :: zalphae(locproma)
    REAL(dp) :: henry,dryreac ! ju_te_20180626

    CHARACTER(len=24) :: trname
    INTEGER  :: jt, jm, naermo, jl, zmode
    INTEGER  :: jtr ! mz_sg_20100320
    LOGICAL  :: l_aerosol(locproma)
    !
    ! check with the value zepdu2 in vdiff
    REAL(dp), PARAMETER :: pepdu2=1.e-3_dp 
    ! smallest radius for which dry deposition velocity is calculated
    REAL(dp), PARAMETER :: crmin=0.01E-8_dp 

    !--- 0) Initialisations: -------------------------------------------------

    zdz(:) = 0._dp
    zzvdrydep(:) = 0._dp

    IF (flag ==I_GP) ti => ti_gp
    IF (flag ==I_LG) ti => ti_lg
    
    ! GENERAL.... BOTH LG AND GP TRACER INDEPENDENT

    ! calc absolute value of wind vector and built max of lower limit
    CALL drydep_calc_landtypefractions(zfrw(1:locproma), &
         zfri(1:locproma), loland(1:locproma), slf(1:locproma,locrow), &
         seaice(1:locproma,locrow))
    
    CALL drydep_calc_vegfrac(zvgrat(1:locproma),zcvbs(1:locproma),  &
         vgrat(1:locproma,locrow),cvs(1:locproma,locrow),           &
         cvw(1:locproma,locrow), loland(1:locproma))
    
    zdz(1:locproma) = deltaZ(_RI_XYZ__(1:locproma,locrow,nlev)) 
    
    abswind(1:locproma) = drydep_calc_abswind(pepdu2   &
         , um1(_RI_XYZ__(1:locproma,locrow,nlev))      &
         , vm1(_RI_XYZ__(1:locproma,locrow,nlev)))
    
    !--- Calculate the aerodynamic resistance:
    CALL drydep_calc_ra(locproma, loland(1:locproma),                  &
         zrahl(1:locproma), zrahw(1:locproma), zrahi(1:locproma),      &
         zrahveg(1:locproma), zrahslsn(1:locproma),                    &
         zustarl(1:locproma), zustarw(1:locproma), zustari(1:locproma),&
         zustveg(1:locproma), zustslsn(1:locproma),                    &
         zdz(1:locproma), cfml(1:locproma,locrow),                     &
         cfncl(1:locproma,locrow),                                     &
         cdnl(1:locproma,locrow), abswind(1:locproma),                 &
         tvir(1:locproma,locrow), tvl(1:locproma,locrow),              &
         ril(1:locproma,locrow), az0(1:locproma,locrow),               &
         z0m(1:locproma,locrow),cdnw(1:locproma,locrow),               &
         cfmw(1:locproma,locrow),cfncw(1:locproma,locrow),             &
         riw(1:locproma,locrow), tvw(1:locproma,locrow),               &
         cdni(1:locproma,locrow), cfmi(1:locproma,locrow),             &
         cfnci(1:locproma,locrow),                                     &
         rii(1:locproma,locrow), tvi(1:locproma,locrow))
       
    !--- Calculate the dry deposition velocity for gases:
       
    !     calling of subroutine in which the "big leaf" dry
    !     deposition velocities are being calculated
    
    abswind(1:locproma) = drydep_calc_abswind(0.01_dp,&
         um1(_RI_XYZ__(1:locproma,locrow,nlev)),vm1(_RI_XYZ__(1:locproma,locrow,nlev)))
    
    CALL drydep_vdbl_parameter(locproma,                  &
         zustveg(1:locproma), tslm1(1:locproma,locrow),   &
         rh_2m(1:locproma,locrow),                        &
         lai(1:locproma,locrow), hc(1:locproma,locrow),   &
         soilph(_RI_XYZ__(1:locproma,locrow,:)))

    ! calculation tracer-dependent 
    ! This subroutine only considers the gas-phase
    ! dry deposition. Aerosol dry deposition is calculated in a
    ! separate routine. The calculations are also only performed
    ! whenever the resistance 'diff' has been defined (that corrects
    ! for differences in diffusivity between the gas and H2O)

    ! ju_te_20180626+
    ! switch on a specific calculation of the cuticular resistance
    ! depending on RH, LAI.. with .NOT.l_ganzeori    
    henry     = 0.0_dp
    dryreac   = 0.0_dp
    ! ju_te_20180626-

    DO jt = 1, DSET(flag)%ntrac
       ! use name of associated regular tracer for tagged tracers
       jtr = ti(jt)%tp%meta%cask_i(I_TAG_REG_IDT)
       trname    = TRIM(ti(jtr)%tp%ident%basename)
       henry     = ti(jt)%tp%meta%cask_r(R_pss)
       dryreac   = ti(jt)%tp%meta%cask_r(R_dryreac_sf)
       IF ((DSET(flag)%idt_vdbl(jt)/= 0) .AND. (DSET(flag)%diff(jt)>0.))  THEN
         IF (.NOT.l_ganzeori) THEN
            CALL drydep_calc_rcut(locproma, rh_2m(1:locproma,locrow),    &
                 zustveg(1:locproma), lai(1:locproma,locrow), dryreac,   &
                 henry,                                                  &
                 rrain_ls(1:locproma,locrow),           & ! ju_te_20190416
                 rrain_cv(1:locproma,locrow),           & ! ju_te_20190416
                 DSET(flag)%rcut_2d(jt)%ptr(1:locproma,locrow),          &
                 DSET(flag)%rws_2d(jt)%ptr(1:locproma,locrow))
            CALL drydep_vdbl(locproma, DSET(flag)%ntrac ,                &
                 loland(1:locproma),                                     &
                 DSET(flag)%zvdrydep(1:locproma,locrow,jt),              &
                 trname, jt,                                             &
                 zfri(1:locproma),    cvs(1:locproma,locrow),            &
                 cvw(1:locproma, locrow), zvgrat(1:locproma),            &
                 zrahw(1:locproma),   zrahi(1:locproma),                 &
                 zrahveg(1:locproma), zrahslsn(1:locproma),              &
                 zustarw(1:locproma), zustari(1:locproma),               &
                 zustveg(1:locproma), zustslsn(1:locproma),              &
                 lai(1:locproma,locrow),                                 &
                 rco_leaf(1:locproma,locrow), fws(1:locproma,locrow),    &
                 DSET(flag)%diff, DSET(flag)%diffrb, DSET(flag)%rmes,    &
                 DSET(flag)%rsoil, DSET(flag)%rwater, DSET(flag)%rsnow,  &
                 rcut_2d=DSET(flag)%rcut_2d(jt)%ptr(1:locproma,locrow),  &
                 rws_2d=DSET(flag)%rws_2d(jt)%ptr(1:locproma,locrow))
         ELSE 
            CALL drydep_vdbl(locproma, DSET(flag)%ntrac,                  &
                 loland(1:locproma),                                      &
                 DSET(flag)%zvdrydep(1:locproma,locrow,jt),               &
                 trname, jt,                                              &
                 zfri(1:locproma),    cvs(1:locproma,locrow),             &
                 cvw(1:locproma, locrow), zvgrat(1:locproma),             &
                 zrahw(1:locproma),   zrahi(1:locproma),                  &
                 zrahveg(1:locproma), zrahslsn(1:locproma),               &
                 zustarw(1:locproma), zustari(1:locproma),                &
                 zustveg(1:locproma), zustslsn(1:locproma),               &
                 lai(1:locproma,locrow),                                  &
                 rco_leaf(1:locproma,locrow), fws(1:locproma,locrow),     &  
                 DSET(flag)%diff, DSET(flag)%diffrb, DSET(flag)%rmes,     &
                 DSET(flag)%rsoil, DSET(flag)%rwater, DSET(flag)%rsnow,   &
                 DSET(flag)%rws, DSET(flag)%rcut)
         END IF
       END IF
    END DO
          

    ! This subroutine only considers the gas-phase
    ! dry deposition. Aerosol dry deposition is calculated in a
    ! separate routine. The calculations are also only performed
    ! whenever the resistance 'diff' has been defined (that corrects
    ! for differences in diffusivity between the gas and H2O)
    
    CALL drydep_vdbl_dealloc
    
    ! start coupling with airsea for land sea mask!
    IF (DSET(flag)%l_airsea_coupling) THEN 
       DO jt = 1, DSET(flag)%ntrac
          IF(DSET(flag)%l_cpl_airsea(jt)) THEN
             DO jl=1,locproma
                ! pure water or ice??
                ! as in airsea (consistent land sea mask!) 
                ! water= (1-land-lake)*(1-seaice)

                IF ((1._dp-slf(jl,locrow)-alake(jl,locrow))&
                     *(1_dp-seaice(jl,locrow)) >= 0.5_dp) THEN
                   DSET(flag)%zvdrydep(jl,locrow,jt)=0.
                  DSET(flag)%vd_bigl(DSET(flag)%idt_vdbl(jt))%ptr(jl,locrow)=0.
                END IF
             END DO
          END IF
       END DO
    END IF
    ! end coupling with airsea for land sea mask!
    
    
    ! CALCULATE AEROSOL DRY DEPOSITION
    if_vdaer: IF (DSET(flag)%l_vdaer) THEN
       
       abswind(1:locproma) = drydep_calc_abswind(0.001_dp,        &
            um1(_RI_XYZ__(1:locproma,locrow,nlev)),               &
            vm1(_RI_XYZ__(1:locproma,locrow,nlev)))
       
       abswind10(1:locproma) = drydep_calc_abswind(0.001_dp,      &
            u10(1:locproma,locrow),v10(1:locproma,locrow)) 
       
       CALL drydep_vdaer_parameter(locproma,            &
            loland(1:locproma),  zustarw(1:locproma),   &
            abswind(1:locproma), abswind10(1:locproma), &
            rh_2m(1:locproma,    locrow),               &
            zalpha(1:locproma),  zalphae(1:locproma),   &
            zbubble(1:locproma), zbeta(1:locproma))
          
       aermo_loop: DO naermo = 1, DSET(flag)%aeromodnum
          
          mod_loop_mass: DO jm =1,DSET(flag)%nmod(naermo)
             
             IF (DSET(flag)%aermodmethod(naermo) == BIN )  CYCLE
             
             radius_mass(:) = 0.
             
             IF (DSET(flag)%l_anymode(naermo)%lptr(jm)) THEN
                l_aerosol(:) =.FALSE.
                DO jl = 1, locproma
                  IF ((DSET(flag)%radii(naermo)%ptr(_RI_XYZN_(jl,locrow,nlev,jm))>crmin)&
                   .AND. (DSET(flag)%densaer(naermo)%ptr(_RI_XYZN_(jl,locrow,nlev,jm)) >&
                        TINY(crmin)) )   l_aerosol(jl) = .TRUE.
                ENDDO
                WHERE (l_aerosol(1:locproma))                     
                   radius_mass(1:locproma) = drydep_calc_radius_mass(  &
                     DSET(flag)%radii(naermo)%ptr(_RI_XYZN_(1:locproma,locrow,nlev,jm)),&
                        DSET(flag)%sigma(naermo)%ptr(jm) )
                ELSEWHERE
                   radius_mass(1:locproma) =  0._dp
                ENDWHERE
                
                ! calculate dry deposition velocities for number
                IF (DSET(flag)%l_anymass(naermo)%lptr(jm))               &
                     CALL drydep_vdaer(zzvdrydep(1:locproma), locproma,    &
                     zvgrat(1:locproma),    cvs(1:locproma, locrow),       &
                     cvw(1:locproma, locrow), zfri(1:locproma),            &
                     zcvbs(1:locproma),     zfrw(1:locproma),              & 
                     zustarw(1:locproma),   zustveg(1:locproma),           & 
                     zustslsn(1:locproma),  zrahl(1:locproma),             &
                     zrahveg(1:locproma),   zrahslsn(1:locproma),          &
                     zevap(1:locproma),     az0(1:locproma, locrow),       &
                     abswind(1:locproma),   radius_mass(1:locproma),       &
                     DSET(flag)%densaer(naermo)%ptr(_RI_XYZN_(1:locproma,locrow,nlev,jm)),&
                     tslm1(1:locproma,locrow),                             &
                     zalpha(1:locproma),    zalphae(1:locproma),           &
                     zbubble(1:locproma),   zbeta(1:locproma))
                
              DSET(flag)%vdrydep_mass(naermo)%ptr(_RI_XYZN_(1:locproma,locrow,1,jm),1)= &
                     zzvdrydep(1:locproma)*1.e-2
                
             ENDIF ! any jm
             
          END DO mod_loop_mass
          
          mod_loop_num: DO jm=1,DSET(flag)%nmod(naermo)
             
             DO jl = 1, locproma
                l_aerosol(jl) =.FALSE.
                IF ((DSET(flag)%radii(naermo)%ptr(_RI_XYZN_(jl,locrow,nlev,jm)) > crmin)&
                 .AND.(DSET(flag)%densaer(naermo)%ptr(_RI_XYZN_(jl,locrow,nlev,jm)) >   &
                     TINY(crmin)) )   l_aerosol(jl) = .TRUE.
             ENDDO
             
             WHERE (l_aerosol(1:locproma))                     
                radius_number(1:locproma) = drydep_calc_radius_num(  &
                     DSET(flag)%radii(naermo)%ptr(_RI_XYZN_(1:locproma,locrow,nlev,jm)))
             ELSEWHERE
                radius_number(1:locproma) =  0._dp
             ENDWHERE
             
             IF (DSET(flag)%l_anynum(naermo)%lptr(jm)                       &
                  .OR. (DSET(flag)%l_anymass(naermo)%lptr(jm).AND.          &
                  DSET(flag)%aermodmethod(naermo) == BIN))                    &
                  CALL drydep_vdaer(zzvdrydep(1:locproma), locproma,          &
                  zvgrat(1:locproma),    cvs(1:locproma, locrow),             &
                  cvw(1:locproma, locrow), zfri(1:locproma),                  &
                  zcvbs(1:locproma),     zfrw(1:locproma),                    &
                  zustarw(1:locproma),   zustveg(1:locproma),                 &
                  zustslsn(1:locproma),  zrahw(1:locproma),                   &
                  zrahveg(1:locproma),   zrahslsn(1:locproma),                &
                  zevap(1:locproma),     az0(1:locproma, locrow),             &
                  abswind(1:locproma),   radius_number(1:locproma),           &
                  DSET(flag)%densaer(naermo)%ptr(_RI_XYZN_(1:locproma,locrow,nlev,jm)),     &
                  tslm1(1:locproma,locrow),                                   &
                  zalpha(1:locproma), zalphae(1:locproma),                    &
                  zbubble(1:locproma), zbeta(1:locproma))
             
             DSET(flag)%vdrydep_num(naermo)%ptr(_RI_XYZN_(1:locproma,locrow,1,jm),1)= &
                  zzvdrydep(1:locproma) *1.e-2
             
          ENDDO mod_loop_num
          
       ENDDO aermo_loop  ! naermo
    END IF if_vdaer
    
    tracer_loop3: DO jt = 1, DSET(flag)%ntrac
       if_drydep: IF(ti(jt)%tp%meta%cask_i(I_drydep) == ON) THEN
          if_Aerosol: IF (ti(jt)%tp%ident%medium == AEROSOL) THEN
             
             naermo = DSET(flag)%traermod(jt)
             IF (naermo == 0) CYCLE
             
             zmode = ti(jt)%tp%meta%cask_i(I_aerosol_mode)
             
             if_bin: IF (DSET(flag)%aermodmethod(naermo) /= BIN ) THEN
                IF (ti(jt)%tp%ident%quantity == AMOUNTFRACTION) &
                     DSET(flag)%zvdrydep(1:locproma,locrow,jt) = &
                     DSET(flag)%vdrydep_mass(naermo)%ptr(&
                     & _RI_XYZN_(1:locproma, locrow, 1, zmode),1)
                IF (ti(jt)%tp%ident%quantity == NUMBERDENSITY) &
                     DSET(flag)%zvdrydep(1:locproma,locrow,jt) = &
                     DSET(flag)%vdrydep_num(naermo)%ptr(&
                     & _RI_XYZN_(1:locproma, locrow, 1, zmode),1)
             ELSE
                DSET(flag)%zvdrydep(1:locproma,locrow,jt) = &
                     DSET(flag)%vdrydep_num(naermo)%ptr(&
                     & _RI_XYZN_(1:locproma, locrow, 1, zmode),1)
             ENDIF if_bin
          END IF if_AEROSOL  !gas already set

          IF (DSET(flag)%idt_vdbl(jt) /= 0) &
          DSET(flag)%vd_bigl(DSET(flag)%idt_vdbl(jt))%ptr(1:locproma,locrow) = &
               100.*DSET(flag)%zvdrydep(1:locproma,locrow,jt)

          zdz(1:locproma) = deltaz(_RI_XYZ__(1:locproma,locrow,nlev))

          DSET(flag)%vd_eff(1:locproma,locrow,jt) = &
               drydep_calc_vd_eff(DSET(flag)%zvdrydep(1:locproma,locrow,jt), &
               time_step_len, zdz(1:locproma) ) 
          
       ENDIF if_drydep
    ENDDO tracer_loop3
    
 END SUBROUTINE calc_vd_drydep
 !==========================================================================

 !==========================================================================
  SUBROUTINE ddep_predef_vdiff

    ! called from ddep_vdiff
    ! calculates dry deposition from predefined surface sink

    USE messy_main_grid_def_mem_bi, ONLY: nlev, nproma, jrow, kproma
    USE messy_main_data_bi,         ONLY: pressi_3d
#if defined (ECHAM5)
    USE messy_main_data_bi,         ONLY: pxtems
#endif
    USE messy_main_timer,         ONLY: time_step_len, delta_time    
#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, pxtm1 => qxtm1 & 
                                      , ti_gp
#else
    USE messy_main_tracer_mem_bi, ONLY: ti_gp
#endif

    ! MESSy
    USE messy_main_tracer,        ONLY: NUMBERDENSITY, AMOUNTFRACTION 
    USE messy_ddep,               ONLY: drydep_posfinit 

    CHARACTER(LEN=*), PARAMETER :: substr='ddep_predef_vdiff'
    REAL(dp) :: zdp(nproma)
#ifdef MESSYTENDENCY
    REAL(dp) :: xtstart(kproma,nlev)
    REAL(dp) :: zxtte(kproma,nlev)
#endif
    REAL(dp) :: xtp1(kproma)
    REAL(dp) :: conv, conv_sum
    REAL(dp) :: zdrydepflux(nproma)
    INTEGER  :: jt
    INTEGER  :: idt
#if defined(ECHAM5)
    REAL(dp), POINTER :: zxtems(:,:)

    zxtems => pxtems(:,1,:,jrow) 
#endif

    zdp(:)           = 0._dp

    zdp(1:kproma) = pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev+1)) - &
                    pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev))

    DO jt=1, N_predepvel  ! loop over number of surface sinks

       idt = X_predepvel(jt)%trid ! tracer ID

       zdrydepflux(:) = 0._dp

#ifndef MESSYTENDENCY
       xtp1(1:kproma) = pxtm1(_RI_X_ZN_(1:kproma, nlev,idt)) + &
                        pxtte(_RI_X_ZN_(1:kproma, nlev,idt)) * time_step_len
#else
       CALL mtend_get_start_l(idt, v0=xtstart)
       xtp1(1:kproma) = xtstart(1:kproma,nlev)
#endif

#if defined(ECHAM5)
       xtp1(1:kproma) = xtp1(1:kproma) + &
            zxtems(1:kproma,idt)/ (zdp(1:kproma)/g)*time_step_len
#endif
       ! xtp1 is in mol/mol
       ! X_predepvel(jt)%field_predepvel is in molec(deposited substance)/m^2/s 
       ! vdiff expects zdrydepflux in mol/mol * kg(air)/m^2/s   
       zdrydepflux(1:kproma) = xtp1(1:kproma)  &
               * X_predepvel(jt)%field_predepvel(1:kproma,jrow) * M_air  &
               * 1.0e-03_dp / N_A 

       ! ensure positive tracer concentrations by adjusting zdrydepflux
       zdrydepflux(1:kproma)= drydep_posfinit( xtp1(1:kproma) &
                              , zdrydepflux(1:kproma) &
                              , time_step_len &
                              , zdp(1:kproma)/g )

       ! sum up tendencies
#if defined (ECHAM5)
       IF (.NOT. l_tendency ) THEN
          zxtems(1:kproma,idt) = zxtems(1:kproma,idt) - zdrydepflux(1:kproma)
       ELSE
#endif
#ifndef MESSYTENDENCY
          pxtte(_RI_X_ZN_(1:kproma, nlev,idt)) = pxtte(_RI_X_ZN_(1:kproma, nlev,idt))  &
               - zdrydepflux(1:kproma)/ (zdp(1:kproma)/g)
#else
          zxtte(:,:) = 0._dp
          zxtte(1:kproma,nlev) = - zdrydepflux(1:kproma)/ (zdp(1:kproma)/g)
          CALL mtend_add_l(my_handle, idt, px=zxtte)
#endif
#if defined (ECHAM5)
       END IF  ! .not.l_tendency
#endif          

       ! diagnostics
       IF (X_predepvel(jt)%import_predepvel%TRLOGICAL) THEN

          SELECT CASE(ti_gp(X_predepvel(jt)%trid)%tp%ident%quantity)
          
          CASE(AMOUNTFRACTION)
             ! mol/mol(tracer) * kg(air)/m^2/s -> molec(tracer)/m^2/s
             conv = conv_mr
             ! mol/mol(tracer) * kg(air)/m^2/s -> kg(tracer)/m^2
             conv_sum = &
               (ti_gp(X_predepvel(jt)%trid)%tp%meta%cask_r(R_molarmass)/M_air) 

          CASE(NUMBERDENSITY)
             conv = conv_nu
             conv_sum = 1.0E-03/M_air

          CASE DEFAULT
             conv = 0.0_dp
             conv_sum = 0.0_dp

          END SELECT

          ! dry deposition flux: molec(tracer)/m^2/s
          X_predepvel(jt)%diagflux(1:kproma,jrow) = &
               zdrydepflux(1:kproma) * conv 

          ! dry deposition: kg(tracer)/m^2
          X_predepvel(jt)%diagsum(1:kproma,jrow)  = &
               X_predepvel(jt)%diagsum(1:kproma,jrow) &
               + zdrydepflux(1:kproma) * delta_time * conv_sum  
       END IF
    END DO

  END SUBROUTINE ddep_predef_vdiff
  !==========================================================================

  !=========================================================================

  SUBROUTINE ddep_read_nml_cpl(status, iou)
   
    ! read namelist for 'coupling' 
    !
    ! Author: Andrea Pozzer, MPICH, Aug 2003

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_main_tracer_mem_bi, ONLY: NGCELL

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit
    ! switch for skipping calculation of Lagrangian 
    ! rate coefficients..it is local,not broadcasted

    NAMELIST /CPL/ L_GP, L_LG, i_lg_method, l_tendency, l_diagflux &
                 , outflux, imp_lai, imp_hc, imp_drag, imp_soilpH  &
                 , rainrate_ls, rainrate_cv                        &
                 , import_predepvel

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='ddep_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1
    outflux = ''


    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

#ifdef ECHAM5
    IF ((L_LG) .AND. (NGCELL > 0)) THEN
       L_LG = .TRUE.
    ELSE
       IF (L_LG) THEN
         WRITE(*,*) 'L_LG = T in namelist'
         WRITE(*,*) 'However no Lagrangian scheme activated ...'
         WRITE(*,*) ' ... setting L_LG = F'
       END IF
       L_LG = .FALSE.
    ENDIF
#endif

    IF (L_GP) THEN
       WRITE(*,*) 'DEPOSITION IN GRIDPOINT SPACE : ON'
    ELSE
       WRITE(*,*) 'DEPOSITION IN GRIDPOINT SPACE : OFF'
    END IF

    IF (L_LG) THEN
       WRITE(*,*) 'DEPOSITION IN LAGRANGIAN SPACE : ON'
       WRITE(*,*) 'DEPOSITION METHOD FOR LG       : ', i_lg_method 
    ELSE
       WRITE(*,*) 'DEPOSITION IN LAGRANGIAN SPACE : OFF'
    END IF

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

    IF (.NOT.L_GP .AND. .NOT.L_LG)  &
       CALL error_bi('  NEITHER GRID POINT NOR LAGRANGIAN REPRESENTATION'//&
       &' REQUESTED (namelist CPL): No calculation possible!' , substr )

  END SUBROUTINE ddep_read_nml_cpl

! ===========================================================================

END MODULE messy_ddep_si
