#include "messy_main_ppd_bi.inc"

! **************************************************************************
MODULE messy_mmd2way_parent_si
! **************************************************************************

  ! MODULE FOR exchanging data with another model of finer resolution
  !
  ! Author: Klaus Ketelsen,  MPICH, Oct 2008
  !         Patrick Joeckel, MPICH, Dec. 2008
  !         Astrid Kerkweg, UNI-MZ, Dec. 2008 -  June 2010
  !                                 Feb. 2012 (unification e5 and c4)
  !                                 2012 - 2016 (extended to 2-Way)
#if defined(MESSYMMD) && !defined(CESM1)
  ! MESSy/BMIL
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                    , debug_bi                         &
                                    , error_bi, info_bi, warning_bi
  ! MESSy/SMCL
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, STRLEN_OBJECT &
                                    , t_chaobj_cpl
  USE messy_main_timer_event,   ONLY: time_event, io_time_event
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, iouerr
  ! MODULE
  USE messy_mmd2way_parent
  USE messy_mmd2way

  IMPLICIT NONE

  INTRINSIC :: NULL
#endif

  PRIVATE
  SAVE

#if defined(MESSYMMD) && !defined(CESM1)
  ! NUMBER OF Childs of this Parent
  INTEGER  :: NumberOfChildren

  LOGICAL, PUBLIC :: lcpl_global_start
  ! -----------------------------------------------------------------------
  !
  ! CHILD CPLDATA STRUCT
  TYPE T_COUPLE_C_DATA
     ! CHANNEL AND CHANNEL OBJECT NAMES IN PARENT
     TYPE(t_chaobj_cpl)                    :: name
     ! POINTER TO PARENT  FIELD
     REAL(DP), POINTER, DIMENSION(:,:,:,:) :: ptr => NULL()
     ! REPRESENTATION OF PARENT FIELD
     INTEGER                               :: rank
     INTEGER, DIMENSION(4)                 :: ldimlen=0
     ! STRING ORDER OF AXES
     CHARACTER(LEN=4)                      :: AXIS
  END TYPE T_COUPLE_C_DATA

  TYPE T_CHILD_DATA
     ! NUMBER OF EXCHANGED FIELDS
     INTEGER                                        :: NEXCH
     TYPE (T_COUPLE_C_DATA), DIMENSION(:), POINTER  :: CPLDATA => NULL()
     ! POINTER for DATA EXCHANGE TIMING
     REAL(DP), POINTER                              :: Waittime => NULL()
     ! LOGICAL for COUPLING TIME STEP (evaluated in init_loop, required in
     !  init_loop and global start)
     LOGICAL                                        :: lcpl = .FALSE.
  END TYPE T_CHILD_DATA

  TYPE (T_CHILD_DATA), DIMENSION(:), ALLOCATABLE    :: CL
  !
  ! -----------------------------------------------------------------------
  ! -----------------------------------------------------------------------
  !
  ! PARENT CPLDATA STRUCT
  TYPE T_EXCH_P_IO
     ! CHANNEL AND CHANNEL OBJECT NAMES IN PARENT
     TYPE(t_chaobj_cpl)                    :: Parentname
     TYPE(t_chaobj_cpl)                    :: Parenttend
     TYPE(t_chaobj_cpl)                    :: Childname
     ! Representation String
     CHARACTER(LEN=STRLEN_MEDIUM)          :: REPR = ''  ! REPRESENTATION STRING
     ! FLAG FOR INTERPOLATION METHOD
     INTEGER                               :: interp = 0 ! INTERPOLATION METHOD
     ! FLAG FOR APPLICATION METHOD overwrite (input field) in GridPoint (appl=0)
     ! FLAG FOR APPLICATION METHOD weighted in GridPoint (1)
     INTEGER                               :: appl = 0   ! application method
     !  weighing factor for "nudging" 1 = hard nudging
     REAL(dp)                              :: fac = 1._dp
  END TYPE T_EXCH_P_IO

  ! MAXIMAL NUMBER OF EXCHANGE FIELDS
  INTEGER, PARAMETER                        :: NMAX_P_EXCH = 100
  TYPE(T_EXCH_P_IO), DIMENSION(NMAX_P_EXCH) :: PFIELD
  TYPE(T_EXCH_P_IO), DIMENSION(NMAX_P_EXCH) :: DEFAULT_PFIELD

  TYPE T_P_EXCHG_DATA_IO
     ! use generalised humidity couplings
     LOGICAL                                :: lgrhum      = .FALSE.
     INTEGER                                :: i_rmy_px    = 0
     REAL(dp)                               :: pcontrol_fi = 30000.
     INTEGER                                :: itype_fw    = 2
     INTEGER                                :: icosexp     = 14
     REAL(dp)                               :: damprel     = 0.2_dp
     INTEGER                                :: itype_VI    = 1
     INTEGER                                :: RCF         = 10000
     INTEGER                                :: RCF_in      = 10000
     LOGICAL                                :: ldiagonly = .FALSE.
     TYPE(T_EXCH_P_IO), DIMENSION(NMAX_P_EXCH) :: FIELD
     ! allow for free boundary layer
     LOGICAL                                :: lfreeslip   = .FALSE.
     LOGICAL                                :: lcpl_gs     = .FALSE.
  END TYPE T_P_EXCHG_DATA_IO

  TYPE(T_P_EXCHG_DATA_IO), DIMENSION(:), ALLOCATABLE :: PIO

  ! PARENT CPLDATA STRUCT
  TYPE T_COUPLE_P_DATA
     ! CHANNEL AND CHANNEL OBJECT NAMES IN PARENT
     TYPE(t_chaobj_cpl)                    :: Parentname
     TYPE(t_chaobj_cpl)                    :: Parenttend
     TYPE(t_chaobj_cpl)                    :: Childname
     ! POINTER TO FIELD in Parent model which will be changed
     ! (i.e. for prognostic variables the tendency !)
     REAL(DP), POINTER, DIMENSION(:,:,:,:) :: target => NULL()
     ! POINTER TO FIELD in Parent model required for tendency calulation
     ! (i.e. the m1 field for prog. vars, otherwise  nothing  !)
     REAL(DP), POINTER, DIMENSION(:,:,:,:) :: targetm1 => NULL()

     ! POINTER TO FIELD yielding the interpolated Child fields
     REAL(DP), POINTER, DIMENSION(:,:,:,:) :: ptr => NULL()
     ! TRACER INDEX
     INTEGER                               :: idt = -99
     ! REPRESENTATION OF PARENT FIELD
     INTEGER                               :: rank
     ! coupling to tendency
     LOGICAL                               :: lte = .FALSE.
     INTEGER, DIMENSION(4)                 :: ldimlen=0
     INTEGER, DIMENSION(4)                 :: gdimlen=0
     ! STRING ORDER OF AXES
     CHARACTER(LEN=4)                      :: AXIS
     ! FLAG FOR INTERPOLATION METHOD
     ! 1 conservative
     ! 2 bilinear ... (to be implemented)

     INTEGER                               :: interp = 0
     ! FLAG FOR APPLICATION METHOD
     !  0: as 1: only possible for input fields (i.e., channel objects made by
     !           'mmd2way_parent' itself, not for tendency calculation not for
     !           tracer  (on the parent side.) the source can of course be a
     !           tracer!)
     !          for appl = 0 the input field is simply transformed to the
     !          full 2-D/ 3-D parent grid (where frac >0.9999) without
     !          application of any weight function
     !     flagged (i.e., they are /= 0 only where the wf is 1
     !  1: field exists only in GridPoint (1) => change tendency
     !     will be worked on in global start
     !  2: dynamical field, treat in init_loop
     !     u,v create => divergence (d), vorticity (vo)
     !     and calculate gradients
     !     ! qqq how to treat temperature (horizontal gradients also required) ?
     INTEGER                               :: appl = 0
     !  weighing factor for "nudging" 1 = hard nudging
     REAL(dp)                              :: fac = 1._dp
     ! Representation String
     CHARACTER(LEN=STRLEN_MEDIUM)          :: REPR
     ! is time dependent ? (relevant Only for COSMO parent)
     LOGICAL                               :: ltimedep
     ! mmd parent has to define its own memory
     LOGICAL                               :: l_sendunit = .FALSE.
  END TYPE T_COUPLE_P_DATA

  TYPE T_COUPLE_P
     ! NUMBER OF EXCHANGE FIELDS
     INTEGER                                       :: NEXCH
     TYPE (T_COUPLE_P_DATA), DIMENSION(:), POINTER :: CPLDATA => NULL()
     ! POINTER for DATA EXCHANGE TIMING
     REAL(DP), POINTER                             :: Waittime => NULL()
     ! Lon /Lat of incoming data
     REAL(kind=dp),      DIMENSION(:,:,:), POINTER :: all_plon => NULL()
     REAL(kind=dp),      DIMENSION(:,:,:), POINTER :: all_plat => NULL()

     ! weight fractions summed over all child PEs
     REAL(dp),             DIMENSION(:,:), POINTER :: frac => NULL()
     REAL(dp),             DIMENSION(:,:), POINTER :: mask => NULL()
     INTEGER                                       :: kmin
     REAL(dp), DIMENSION(:), POINTER               :: kminfac => NULL()

     ! weighting function
     REAL(dp),             DIMENSION(:,:), POINTER :: wf => NULL()

     ! ENTRIES FOR GENERALIZED HUMIDITY COUPLING:
     ! use generalised coupling
     LOGICAL                                       :: lgrhum = .FALSE.
     ! data index for vapor
     INTEGER                                       :: idx_qv = -99
     ! data index for cloud water
     INTEGER                                       :: idx_qc = -99
     ! data index for cloud ice
     INTEGER                                       :: idx_qi = -99
     ! data index for meridional wind
     INTEGER                                       :: idx_u  = -99
     ! data index for zonal wind
     INTEGER                                       :: idx_v  = -99
     ! data index for deriv. of meridional wind
     INTEGER                                       :: idx_t    = -99
     ! data index for surface pressure
     INTEGER                                       :: idx_ps   = -99

     ! NAMELIST SWITCHES

     ! treatment of boundary zone for back transition
     INTEGER                                       :: i_rmy_px  = 0
     REAL(dp)                                      :: pcontrol_fi = 30000.

     INTEGER                                       :: itype_fw  = 2
     INTEGER                                       :: icosexp   = 14
     REAL(dp)                                      :: damprel   = 0.2_dp
     INTEGER                                       :: RCF       = 10000
     INTEGER                                       :: RCF_in    = 10000
     LOGICAL                                       :: ldiagonly = .FALSE.
     ! itype_VI switches the vertical interpolation  of C2E
     ! == 1  => NCREGRID
     ! == 2  => INT2LM inverse (MesoTel version)
     INTEGER                                       :: itype_VI  = 1
     !
     ! allow for not forced boundary layer
     LOGICAL                                       :: lfreeslip = .FALSE.
     LOGICAL                                       :: lcpl_gs   = .FALSE.

     REAL(DP), POINTER, DIMENSION(:,:) :: aps => NULL()

  END type T_COUPLE_P

  TYPE(T_COUPLE_P), DIMENSION(:), ALLOCATABLE :: PAR

  ! -----------------------------------------------------------------------

  ! TIMER FOR COUPLING
  TYPE(time_event), DIMENSION(:), ALLOCATABLE    :: CPL_EVENT
  TYPE(io_time_event), DIMENSION(:), ALLOCATABLE :: CPL_IOEVENT

  ! -----------------------------------------------------------------------
  LOGICAL, PARAMETER :: lparent_numbering = .FALSE.
  LOGICAL, PARAMETER :: ldev              = .FALSE.
  !LOGICAL, PARAMETER :: ldev              = .TRUE.

  ! INTERFACES
  INTERFACE mmd2way_parent_initialize
    MODULE PROCEDURE mmd2way_parent_initialize
  END INTERFACE mmd2way_parent_initialize

  INTERFACE mmd2way_parent_init_memory
    MODULE PROCEDURE mmd2way_parent_init_memory
  END INTERFACE mmd2way_parent_init_memory

  INTERFACE mmd2way_parent_init_coupling
    MODULE PROCEDURE mmd2way_parent_init_coupling
  END INTERFACE mmd2way_parent_init_coupling

  INTERFACE mmd2way_parent_global_start
    MODULE PROCEDURE mmd2way_parent_global_start
  END INTERFACE mmd2way_parent_global_start

  INTERFACE mmd2way_parent_global_end
    MODULE PROCEDURE mmd2way_parent_global_end
  END INTERFACE mmd2way_parent_global_end

  INTERFACE mmd2way_parent_free_memory
    MODULE PROCEDURE mmd2way_parent_free_memory
  END INTERFACE mmd2way_parent_free_memory

  PUBLIC :: mmd2way_parent_initialize
  PUBLIC :: mmd2way_parent_init_coupling
  PUBLIC :: mmd2way_parent_init_memory

  PUBLIC :: mmd2way_parent_global_start
  PUBLIC :: mmd2way_parent_global_end
  PUBLIC :: mmd2way_parent_free_memory

CONTAINS

! ---------------------------------------------------------------------
  SUBROUTINE mmd2way_parent_initialize

    ! MESSy/BMIL
    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    ! MMD
    USE mmd_parent,               ONLY: MMD_P_Init            &
                                      , MMD_Parent_for_Child  &
                                      , MMD_P_Allocate_Child
    ! MESSy/SMCL
    USE messy_main_tools,         ONLY: find_next_free_unit, int2str

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! LOCAL
    INTEGER                     :: ic             ! loop over childs
    INTEGER                     :: child_id       ! child index
    CHARACTER(LEN=*), PARAMETER :: substr = 'mmd2way_parent_initialize'
    CHARACTER(LEN=3)            :: NML_TAG
    INTEGER                     :: iou            ! unit for namelist file
    INTEGER                     :: status         ! error status
    LOGICAL                     :: lgrh            = .FALSE.
    LOGICAL                     :: lgrh_def        = .FALSE.
    INTEGER                     :: i_rmy_def       = 0
    INTEGER                     :: i_rmy_px        = 0
    REAL(dp)                    :: pcontrol_fi     = 30000._dp
    REAL(dp)                    :: pcon_fi_def     = 30000._dp
    LOGICAL                     :: l_freeslip      = .FALSE.
    LOGICAL                     :: l_freeslip_def  = .FALSE.
    LOGICAL                     :: lcpl_gs         = .FALSE.
    LOGICAL                     :: lcpl_gs_def     = .FALSE.
    INTEGER                     :: itype_fw        = 2
    INTEGER                     :: itype_fw_def    = 2
    INTEGER                     :: icosexp         = 14
    INTEGER                     :: icosexp_def     = 14
    REAL(dp)                    :: damprel         = 0.2_dp
    REAL(dp)                    :: damprel_def     = 0.2_dp
    INTEGER                     :: itype_VI        = 1
    INTEGER                     :: itype_VI_def    = 1
    INTEGER                     :: RCF             = 10000
    INTEGER                     :: RCF_in          = 10000
    INTEGER                     :: RCF_def         = 10000
    INTEGER                     :: RCF_in_def      = 10000
    LOGICAL                     :: ldiagonly       = .FALSE.
    LOGICAL                     :: ldiagonly_def   = .FALSE.

    CALL start_message_bi(submodstr, 'INITIALIZE MMD', substr)

    NumberOfChildren = SIZE(MMD_Parent_for_Child)

    IF (ldev) THEN
       CALL MMD_P_Allocate_Child(NumberOfChildren, l2way=.TRUE.)
    ELSE
       CALL MMD_P_Allocate_Child(NumberOfChildren, l2way=.FALSE.)
    ENDIF
    DO ic=1,NumberOfChildren
       child_id = MMD_Parent_for_Child(ic)
       CALL MMD_P_Init (child_id, ic)
    END DO

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! PARENT COUPLING
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! DIMENSION STRUCT containing the namelist definition for CPLDATA
    ! TO the NUMBER OF CHILDS (CHILDs)
    IF (ldev) ALLOCATE(PAR(NumberOfChildren))
    ALLOCATE(PIO(NumberOfChildren))

    ! INITIALISE NAMELIST SETTINGS

    IF (ldev) THEN
    ! 1.) read default namelist
       CALL warning_bi('First namelist is used as default for coupling', substr)
       IF (p_parallel_io) THEN
          iou = find_next_free_unit(100,200)

       CALL mmd2way_parent_read_nml_cpl_child(status, iou, 'DEF'              &
            , lgrh_def, i_rmy_def                                             &
            , pcon_fi_def, l_freeslip_def, lcpl_gs_def, itype_fw_def          &
            , icosexp_def, damprel_def, itype_VI_def, RCF_def, RCF_IN_DEF     &
            , ldiagonly_def)

       IF (status /= 0) CALL error_bi( &
            'error in read CPL_PAR_CHILD namelist ',substr)

       DEFAULT_PFIELD = PFIELD

       lcpl_global_start = lcpl_gs_def

       ! 2.) Child specific settings
       child_loop_01: DO ic=1,NumberOfChildren
          IF (lparent_numbering) THEN
             CALL int2str(NML_TAG, ic)
          ELSE
             CALL int2str(NML_TAG, MMD_Parent_for_Child(ic))
          ENDIF

          CALL mmd2way_parent_read_nml_cpl_child(                             &
                 status, iou, NML_TAG, lgrh, i_rmy_px                         &
               , pcontrol_fi, l_freeslip, lcpl_gs, itype_fw, icosexp, damprel &
               , itype_VI, RCF, RCF_IN, ldiagonly)

          IF (status == 0) THEN
             ! USE INSTANCE SPECIFIC NAMELIST
             CALL info_bi(&
                  'using instance specific namelist for instance '//NML_TAG &
                  , substr)
             PIO(ic)%lgrhum      = lgrh
             PIO(ic)%i_rmy_px    = i_rmy_px
             PIO(ic)%pcontrol_fi = pcontrol_fi
             PIO(ic)%FIELD       = PFIELD
             PIO(ic)%lfreeslip   = l_freeslip
             PIO(ic)%lcpl_gs     = lcpl_gs
             PIO(ic)%itype_fw    = itype_fw
             PIO(ic)%icosexp     = icosexp
             PIO(ic)%damprel     = damprel
             PIO(ic)%itype_VI    = itype_VI
             PIO(ic)%rcf         = rcf
             PIO(ic)%rcf_in      = rcf_in
             PIO(ic)%ldiagonly   = ldiagonly

             lcpl_global_start   =  lcpl_global_start .OR. PIO(ic)%lcpl_gs
          ELSE IF (status == 444) THEN
             ! USE DEFAULT
             CALL info_bi(&
                  'using default namelist settings for instance '//NML_TAG &
                  , substr)
             PIO(ic)%lgrhum      = lgrh_def
             PIO(ic)%i_rmy_px    = i_rmy_def
             PIO(ic)%pcontrol_fi = pcon_fi_def
             PIO(ic)%FIELD       = DEFAULT_PFIELD
             PIO(ic)%lfreeslip   = l_freeslip_def
             PIO(ic)%lcpl_gs     = lcpl_gs_def
             PIO(ic)%itype_fw    = itype_fw_def
             PIO(ic)%icosexp     = icosexp_def
             PIO(ic)%damprel     = damprel_def
             PIO(ic)%itype_VI    = itype_VI_def
             PIO(ic)%rcf         = rcf_def
             PIO(ic)%rcf_in      = rcf_in_def
             PIO(ic)%ldiagonly   = ldiagonly_def
          ELSE
             CALL error_bi(&
                  'error in read CPL_PAR_CHILD namelist '//NML_TAG, substr)
          ENDIF

       END DO child_loop_01
    ENDIF ! p_parallel_io
    ENDIF ! ldev
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! CHILD COUPLING
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! ALLOCATE CPL_EVENT FOR CHILDS
    ALLOCATE(CL(NumberOfChildren))
    ALLOCATE(CPL_IOEVENT(NumberOfChildren))
    ALLOCATE(CPL_EVENT(NumberOfChildren))

    CALL end_message_bi(submodstr, 'INITIALIZE MMD', substr)

    RETURN

  END SUBROUTINE mmd2way_parent_initialize
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE mmd2way_parent_init_memory

    ! MESSy/BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: SCALAR, GP_2D_HORIZONTAL, GP_1D_ILEV
    ! MESSy/SMCL
    USE messy_main_channel,      ONLY: new_channel, new_channel_object &
                                     , new_attribute
    USE messy_main_tools,        ONLY: int2str
    ! MMD
    USE mmd_parent,              ONLY: MMD_Parent_for_Child

    IMPLICIT NONE

    ! LOCAL
    INTEGER                     :: ic
    INTEGER                     :: child_id
    INTEGER                     :: status
    CHARACTER(LEN=2)            :: clstr
    CHARACTER(LEN=*), PARAMETER :: substr = 'mmd2way_parent_init_memory'

    CALL start_message_bi(submodstr, 'INIT MEMORY MMD2WAY_PARENT', substr)

    ! DEFINE SUBMODEL-CHANNEL
    CALL new_channel(status, submodstr)
    CALL channel_halt(substr, status)

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! CHILD COUPLING
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    DO ic=1, NumberOfChildren

       child_id = MMD_Parent_for_Child(ic)
       CALL Setup_Child_Timer (child_id,ic)

       CALL int2str(clstr, child_id)

       ! ALLOCATE TIME ARRAY FOR PARENT GET
       CALL new_channel_object(status, submodstr, 'WaitTime_'//clstr &
            , p0=CL(ic)%WaitTime, reprid = SCALAR)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, submodstr, 'WaitTime_'//clstr &
            , 'long_name', c='time Parent waits for Child'//clstr)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, submodstr, 'WaitTime_'//clstr &
            , 'units', c='s')
       CALL channel_halt(substr, status)

       IF (ldev ) THEN
       ! weight fraction channels: weight fractions summed over all child PEs
       CALL new_channel_object(status, submodstr, 'weightfrac_'//clstr &
            , p2=PAR(ic)%frac, reprid=GP_2D_HORIZONTAL)
       CALL channel_halt(substr, status)

       ! weight fraction channels: weight fractions summed over all child PEs
       CALL new_channel_object(status, submodstr, 'mask_'//clstr &
            , p2=PAR(ic)%mask, reprid=GP_2D_HORIZONTAL)
       CALL channel_halt(substr, status)

       ! horizontal weight function
       CALL new_channel_object(status, submodstr, 'weightfunc_'//clstr &
            , p2=PAR(ic)%wf, reprid=GP_2D_HORIZONTAL)
       CALL channel_halt(substr, status)

       ! vertical weight function
       CALL new_channel_object(status, submodstr, 'kminfac_'//clstr &
            , p1=PAR(ic)%kminfac, reprid=GP_1D_ILEV)
       CALL channel_halt(substr, status)

       ENDIF

     END DO

     CALL end_message_bi(modstr, 'INIT MEMORY MMD2WAY_PARENT', substr)

  END SUBROUTINE mmd2way_parent_init_memory
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE mmd2way_parent_init_coupling

    ! MESSy/BMIL
    USE messy_main_grid_def_mem_bi,  ONLY: nlev, nproma, ngpblks
    USE messy_main_grid_def_bi,      ONLY: hyam, hybm
    USE messy_main_mpi_bi,           ONLY: p_pe, p_io, p_bcast
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL
#ifdef ECHAM5
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_tracer,           ONLY: get_tracer, full2base_sub
#endif
    ! MESSy/SMCL
    USE messy_main_channel_repr, ONLY: get_representation_info &
                                     , get_representation_id
    USE messy_main_channel,      ONLY: get_channel_object      &
                                     , get_channel_object_info &
                                     , new_channel_object      &
                                     , new_attribute
    USE messy_main_tools,        ONLY: int2str

    ! MMD
    USE mmd_parent,              ONLY: MMD_Parent_for_Child          &
                                     , MMD_P_Get_DataArray_Name      &
                                     , MMD_P_GetNextParArray         &
                                     , MMD_P_Set_ParDataArray        &
                                     , MMD_P_setInd_and_AllocMem     &
                                     , MMD_Recv_from_Child           &
                                     , MMD_STATUS_OK                 &
                                     , MMD_P_Set_ParDataArray_Name   &
                                     , MMD_P_Set_ParDataArray_EndList&
                                     , MMD_Inter_Bcast

    ! ********** height reduced representation **********
    ! MESSy/BMIL
    USE messy_main_channel_bi,           ONLY: DIMID_LON, DIMID_LAT      &
                                             , DIMID_TLV                 &
                                             , DC_GP, gp_nseg            &
                                             , gp_start, gp_cnt          &
                                             , gp_meml, gp_memu


    ! MESSy/SMCL
    USE messy_main_channel_dimensions,   ONLY: new_dimension             &
                                             , add_dimension_variable    &
                                             , add_dimension_variable_att
    USE messy_main_channel_repr,         ONLY: new_representation        &
                                             , set_representation_decomp &
                                             , AUTO, IRANK, PIOTYPE_COL  &
                                             , repr_def_axes
    USE messy_main_channel,              ONLY: DIMID_UNDEF, REPR_UNDEF
    USE messy_main_constants_mem,        ONLY: pi, STRLEN_MEDIUM

    ! ********** height reduced representation **********

    IMPLICIT NONE

    INTRINSIC :: INT, TRIM

    ! LOCAL
#ifdef COSMO
    CHARACTER(LEN=*), PARAMETER     :: substr='cosmo_parent_init_coupling'
#else
    CHARACTER(LEN=*), PARAMETER     :: substr='echam_parent_init_coupling'
#endif
    INTEGER                          :: ic
    INTEGER                          :: child_id
    INTEGER                          :: ii
    CHARACTER(LEN=STRLEN_CHANNEL)    :: myChannel
    CHARACTER(LEN=STRLEN_OBJECT)     :: myName
    INTEGER                          :: status, status_tmp
    INTEGER                          :: reprid
    CHARACTER(LEN=2)                 :: clstr = ''
#ifdef ECHAM5
    CHARACTER(LEN=STRLEN_MEDIUM)     :: trbasename, chbasename
    CHARACTER(LEN=STRLEN_MEDIUM)     :: trsubname, chsubname
#endif

    REAL(kind=dp), DIMENSION(1)      :: p_rdheight
    REAL(kind=dp), PARAMETER         :: p0sl=1.e5
    REAL(kind=dp)                    :: p0hl
    INTEGER                          :: k, nk, kered

    INTEGER, DIMENSION(1)            :: ikmin
    CHARACTER(LEN=3)                 :: IC_TAG = ''
    INTEGER, DIMENSION(10)           :: DIMID_KMINX       = DIMID_UNDEF
    INTEGER, DIMENSION(10)           :: GP_3D_RDCX        = REPR_UNDEF
    INTEGER, DIMENSION(10)           :: DIMID_KMINXp1     = DIMID_UNDEF
    INTEGER, DIMENSION(10)           :: GP_3D_RDCXp1      = REPR_UNDEF

    INTEGER, DIMENSION(11)           :: intbuffer
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    INTEGER                          :: dimids(4)
    INTEGER                          :: ix

    REAL(DP), ALLOCATABLE, DIMENSION(:) :: array
    CHARACTER(LEN=STRLEN_MEDIUM)        :: unit = ''

    CALL start_message_bi(modstr, 'INIT Child COUPLING', substr)

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! PARENT COUPLING
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF (ldev) THEN
       ! INTERPRET NAMELISTS FOR COULING WITH CHILDs
       CALL start_message_bi(modstr, 'INTERPRET NAMELIST', substr)
       CALL interpret_parent_namelist
       CALL end_message_bi(modstr, 'INTERPRET NAMELIST', substr)

       ! Set Names of Data Arrays for coupling IN MMD
       status = MMD_STATUS_OK

       DO ic = 1, NumberOfChildren
          DO ii=1, Par(ic)%NEXCH
             IF (Par(ic)%CPLDATA(ii)%Parentname%CHA(1:7) == 'mmd2way') THEN
                IF (Par(ic)%CPLDATA(ii)%Childname%CHA(1:7) /= 'mmd2way') &
                     Par(ic)%CPLDATA(ii)%l_sendunit = .TRUE.
             ELSE
                Par(ic)%CPLDATA(ii)%l_sendunit = .FALSE.
             ENDIF
             CALL MMD_P_Set_ParDataArray_Name( ic &
                  , TRIM(Par(ic)%CPLDATA(ii)%Parentname%CHA) &
                  , TRIM(Par(ic)%CPLDATA(ii)%Parentname%OBJ) &
                  , TRIM(Par(ic)%CPLDATA(ii)%Childname%CHA)  &
                  , TRIM(Par(ic)%CPLDATA(ii)%Childname%OBJ)  &
                  , TRIM(Par(ic)%CPLDATA(ii)%REPR)           &
                  , Par(ic)%CPLDATA(ii)%interp               &
                  , Par(ic)%CPLDATA(ii)%l_sendunit           &
                  , status_tmp)
             IF(status_tmp /= MMD_STATUS_OK) status = status_tmp
          ENDDO

          IF(status /= MMD_STATUS_OK)   THEN
             CALL int2str(clstr, ic)
             CALL error_bi('Error in MMD_Set_ParDataArray_Name for Child'//clstr &
                  ,substr)
          ENDIF

          ! CLOSE LIST
          CALL MMD_P_Set_ParDataArray_EndList(ic)

          ! exchange information about humidity coupling
          child_id = MMD_Parent_for_Child(ic)
          ! use generalised humidity?
          IF (PAR(ic)%lgrhum) THEN
             intbuffer(1) = 1
          ELSE
             intbuffer(1) = 0
          ENDIF
          ! exclude additional grid boxes at domain border
          intbuffer(2) = PAR(ic)%i_rmy_px
          ! pcontrol_fi
          intbuffer(3) = INT(PAR(ic)%pcontrol_fi)
          ! type der Relaxierung
          intbuffer(4) = PAR(ic)%itype_fw
          ! Cosinus exponent for itype_fw
          intbuffer(5) = PAR(ic)%icosexp
          ! relativ length of cosinus function, here
          intbuffer(6) = INT(PAR(ic)%damprel * 100)
          ! type der vertical interpolation for C2E
          ! == 1 => NCREGRID  == 2 INT2LM inverse
          intbuffer(7) = PAR(ic)%itype_VI
          ! perform coupling of Child data to parent in global_start (1)
          ! DEFAULT: global_end (0)
          intbuffer(8) = 0
          IF (lcpl_global_start) intbuffer(8) = 1
          intbuffer(9)  = PAR(ic)%RCF
          intbuffer(10) = PAR(ic)%RCF_IN
          IF (PAR(ic)%ldiagonly) THEN
             intbuffer(11) = 1
          ELSE
             intbuffer(11) = 0
          ENDIF

          CALL MMD_Inter_Bcast (intbuffer, .TRUE.,child_id)
       END DO ! loop over childs
    ENDIF ! ldev
    ! mmd2way_parent MODULE ROUTINE (ECHAM-5 / COSMO INTERFACE)
    !
    DO ic=1,NumberOfChildren
       child_id = MMD_Parent_for_Child(ic)

       IF (ldev) THEN
          IF (p_pe == 0)  THEN
             CALL MMD_Recv_from_Child (child_id, p_rdheight,1, 0, 110, status)
             IF (status /= 0) &
                  CALL error_bi('MMD_Recv_from_Child ERROR 110', substr)

             ! calculate PARENT level index, above which coupling should be zero
             DO k=SIZE(hyam),1,-1
                p0hl = hyam(k) + hybm(k) * p0sl
                IF(p0hl <= p_rdheight(1)) THEN
                   Par(ic)%kmin = k
                   EXIT
                END IF
             END DO
          END IF
          CALL p_bcast(Par(ic)%kmin, p_io)

          ikmin = Par(ic)%kmin
          CALL MMD_Inter_Bcast (ikmin, .TRUE.,child_id)

       ENDIF ! ldev

       CALL MMD_P_Get_DataArray_Name(ic)
       CALL Setup_Child_Area (child_id,ic)
       IF (ldev) CALL Setup_ParData_Exchange(child_id,ic) ! Parent Coupling
       CALL Setup_Data_Exchange_with_Child (child_id,ic)
       CALL Define_data_arrays (child_id,ic)
    END DO

    CALL end_message_bi(modstr, 'INIT Child COUPLING', substr)

    CALL start_message_bi(modstr, 'INIT Parent COUPLING', substr)
    ! ---------------------------------------------------------------------
    ! PARENT COUPLING +
    IF (ldev) THEN

    DO ic=1,NumberOfChildren
       IF (Par(ic)%NEXCH /= 0) THEN
          child_id = MMD_Parent_for_Child(ic)

          ! ---------------------------------------------
          kered = nlev-Par(ic)%kmin+1
          CALL int2str(IC_TAG, ic)
          CALL new_dimension(status, DIMID_KMINX(ic), 'lev_rdc'//IC_TAG      &
               , kered)
          CALL channel_halt(substr, status)
          ALLOCATE(array(kered))
          DO ix=Par(ic)%kmin, nlev
             array(ix-Par(ic)%kmin+1) = REAL(ix, DP)
          END DO
          CALL add_dimension_variable(status, 'lev_rdc'//IC_TAG &
               , 'lev_rdc'//IC_TAG , array)
          CALL channel_halt(substr, status)
          CALL add_dimension_variable_att(status, 'lev_rdc'//IC_TAG &
               , 'lev_rdc'//IC_TAG , 'positive', c='down')
          CALL channel_halt(substr, status)
          DEALLOCATE(array)

          CALL new_representation(status, GP_3D_RDCX(ic), 'GP_3D_RDCX'//IC_TAG &
               , rank = 3, link = 'xxx-', dctype = DC_GP                       &
               , dimension_ids = (/ &
                 _RI_XYZ__(DIMID_LON, DIMID_LAT, DIMID_KMINX(ic)) /)           &
               , ldimlen       = (/ _RI_XYZ__(nproma , ngpblks, AUTO) /)       &
               , output_order  = (/ _IX_XYZ__,_IY_XYZ__,_IZ_XYZ__ /)           &
               , axis = repr_def_axes(_RI_XYZ__('X','Y','Z'),'-')              &
               )
          CALL channel_halt(substr, status)

          nseg = gp_nseg
          ALLOCATE(start(nseg,IRANK))
          ALLOCATE(cnt(nseg,IRANK))
          ALLOCATE(meml(nseg,IRANK))
          ALLOCATE(memu(nseg,IRANK))

          start(:,:) = gp_start(:,:)
          cnt(:,:)   = gp_cnt(:,:)
          meml(:,:)  = gp_meml(:,:)
          memu(:,:)  = gp_memu(:,:)

          cnt(:,_IZ_XYZ__)  = kered
          memu(:,_IZ_XYZ__) = kered

          CALL set_representation_decomp(status, GP_3D_RDCX(ic) &
               , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
          CALL channel_halt(substr, status)

          DEALLOCATE(start) ; NULLIFY(start)
          DEALLOCATE(cnt)   ; NULLIFY(cnt)
          DEALLOCATE(meml)  ; NULLIFY(meml)
          DEALLOCATE(memu)  ; NULLIFY(memu)

          ! ---------------------------------------------

          ! ---------------------------------------------
          kered = nlev-Par(ic)%kmin+1 +1
          CALL int2str(IC_TAG, ic)
          CALL new_dimension(status, DIMID_KMINXp1(ic), 'lev_rdcp1'//IC_TAG &
               , kered)
          CALL channel_halt(substr, status)
          ALLOCATE(array(kered))
          DO ix=Par(ic)%kmin, nlev+1
             array(ix-Par(ic)%kmin+1) = REAL(ix, DP)
          END DO
          CALL add_dimension_variable(status, 'lev_rdcp1'//IC_TAG     &
               , 'lev_rdcp1'//IC_TAG , array)
          CALL channel_halt(substr, status)
          CALL add_dimension_variable_att(status, 'lev_rdcp1'//IC_TAG &
               , 'lev_rdcp1'//IC_TAG , 'positive', c='down')
          CALL channel_halt(substr, status)
          DEALLOCATE(array)

          CALL new_representation(status, GP_3D_RDCXp1(ic)                     &
               , 'GP_3D_RDCXp1'//IC_TAG                                        &
               , rank = 3, link = 'xxx-', dctype = DC_GP                       &
               , dimension_ids = (/ &
                 _RI_XYZ__(DIMID_LON, DIMID_LAT, DIMID_KMINXp1(ic))/)          &
               , ldimlen       = (/ _RI_XYZ__(nproma , ngpblks, AUTO) /)       &
               , output_order  = (/ _IX_XYZ__,_IY_XYZ__,_IZ_XYZ__ /)           &
               , axis = repr_def_axes(_RI_XYZ__('X','Y','Z'),'-')              &
               )
          CALL channel_halt(substr, status)


          nseg = gp_nseg
          ALLOCATE(start(nseg,IRANK))
          ALLOCATE(cnt(nseg,IRANK))
          ALLOCATE(meml(nseg,IRANK))
          ALLOCATE(memu(nseg,IRANK))

          start(:,:) = gp_start(:,:)
          cnt(:,:)   = gp_cnt(:,:)
          meml(:,:)  = gp_meml(:,:)
          memu(:,:)  = gp_memu(:,:)

          cnt(:,_IZ_XYZ__)  = kered
          memu(:,_IZ_XYZ__) = kered

          CALL set_representation_decomp(status, GP_3D_RDCXp1(ic) &
               , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
          CALL channel_halt(substr, status)

          DEALLOCATE(start) ; NULLIFY(start)
          DEALLOCATE(cnt)   ; NULLIFY(cnt)
          DEALLOCATE(meml)  ; NULLIFY(meml)
          DEALLOCATE(memu)  ; NULLIFY(memu)

          ! ---------------------------------------------
          ! decrease coupling at (CHILD) top with cosine function

          Par(ic)%kminfac(:)               = 1._dp
          Par(ic)%kminfac(1:Par(ic)%kmin)  = 0.0_dp
          !
          nk = 10 ! smoothing over nk levels
          DO k = 1,nk
             Par(ic)%kminfac(Par(ic)%kmin+k) = &
                  COS((pi/2._dp)-(pi/2._dp)/(nk-1)*(k-1))**2._dp
          END DO

          IF (PAR(ic)%lfreeslip) THEN
             ! let lowest 3 level free
             Par(ic)%kminfac(nlev-2:nlev+1) = 0._dp

             nk = 3
             DO k = 0, nk-1
                Par(ic)%kminfac(nlev-k) = &
                     COS((pi/2._dp)-(pi/2._dp)/(nk-1)*k)**2._dp
             END DO
          END IF

          ii = 0
          DO WHILE (MMD_P_GetNextParArray (ic, myChannel, myName))
             ii = ii + 1

                ! Make memory for data exchange field
                ! 1. get representation id. Use namelist entry
                CALL get_representation_id(status &
                     , Par(ic)%CPLData(ii)%repr, reprid)
                CALL channel_halt (substr, status)

#ifdef ECHAM5
                ! 2. find target field / make target memory
                if_tracer: IF (TRIM(MyChannel) == 'tracer_gp' .OR. &
                     TRIM(MyChannel) == 'tracer_gp_te') THEN
                   Par(ic)%CPLDATA(ii)%lte =.TRUE.
                   CALL full2base_sub(status &
                        , TRIM(Par(ic)%CPLDATA(ii)%Parentname%CHA)&
                        , chbasename, chsubname)
                   CALL full2base_sub(status &
                        , TRIM(Par(ic)%CPLDATA(ii)%Parentname%OBJ)&
                        , trbasename, trsubname)
                   CALL get_tracer( status, TRIM(chsubname), TRIM(trbasename) &
                        , TRIM(trsubname), idx=Par(ic)%CPLData(ii)%idt)
                   CALL tracer_halt(substr,status)
                   ! make memory
                   CALL int2str(clstr, child_id)
                   CALL new_channel_object(status, submodstr        &
                        , TRIM(myName)//'_t'//clstr                 &
                        , p4=Par(ic)%CPLDATA(ii)%target             &
                        , reprid = reprid, lrestreq=.TRUE.)
                   CALL channel_halt (substr, status)
                   Par(ic)%CPLDATA(ii)%target = 1.e-35
                ELSE ! if_tracer
#endif
                   if_mmd2way: IF (MyChannel(1:7) /= 'mmd2way') THEN
                      if_lte: IF (TRIM(Par(ic)%CPLDATA(ii)%Parenttend%OBJ) &
                           /= '') THEN

                         Par(ic)%CPLData(ii)%lte =.TRUE.

                         CALL get_channel_object (status &
                              , TRIM(Par(ic)%CPLDATA(ii)%Parenttend%CHA) &
                              , TRIM(Par(ic)%CPLDATA(ii)%Parenttend%OBJ) &
                              , p4=Par(ic)%CPLData(ii)%target)
                         IF (status /=0)  write(iouerr,*)                &
                          'MMD2WAY_PARENT ERROR: Target tendency not available'&
                              , TRIM(MyChannel),' ',TRIM(myName)
                         CALL channel_halt(substr, status)
                         CALL get_channel_object (status     &
                              , TRIM(MyChannel), TRIM(myName)&
                              , p4=Par(ic)%CPLData(ii)%targetm1)
                         IF (status /=0)  write(iouerr,*)                 &
                              'MMD2WAY_PARENT ERROR: Target not available'&
                              , TRIM(MyChannel),' ',TRIM(myName)
                         CALL channel_halt(substr, status)

                         ! b) get dimension id info
                         CALL get_channel_object_info(status                  &
                              , TRIM(MyChannel), TRIM(myName), reprid = reprid)
                         CALL channel_halt(substr, status)
                         CALL get_representation_info(status,' ' &
                              , ID=reprid, dim_ids=dimids)
                         CALL channel_halt(substr, status)
                         Par(ic)%CPLData(ii)%ltimedep = .FALSE.
                         DO ix=1,4
                            IF (dimids(ix) ==  DIMID_TLV) THEN
                               Par(ic)%CPLData(ii)%ltimedep = .TRUE.
                               EXIT
                            ENDIF
                         END DO

                      ELSE ! lte
                         CALL get_channel_object (status      &
                              , TRIM(MyChannel), TRIM(myName) &
                              , p4=Par(ic)%CPLData(ii)%target )
                         IF (status /=0)  write(iouerr,*) &
                              'MMD2WAY_PARENT ERROR: Target not available'&
                              , TRIM(MyChannel),' ',TRIM(myName)
                         CALL channel_halt(substr, status)
                         NULLIFY(Par(ic)%CPLData(ii)%targetm1)

                         ! b) get dimension id info
                         CALL get_channel_object_info(status  &
                              , TRIM(MyChannel), TRIM(myName) &
                              ,reprid = reprid)
                         CALL channel_halt(substr, status)
                         CALL get_representation_info(status,' ' &
                              , ID=reprid, dim_ids=dimids)
                         CALL channel_halt(substr, status)

                         Par(ic)%CPLData(ii)%ltimedep = .FALSE.
                         DO ix=1,4
                            IF (dimids(ix) ==  DIMID_TLV) THEN
                               Par(ic)%CPLData(ii)%ltimedep = .TRUE.
                               EXIT
                            ENDIF
                         END DO
                      END IF if_lte
                   ELSE ! mmd2way
                      ! make memory
                      CALL int2str(clstr, child_id)
                      CALL new_channel_object(status, submodstr        &
                           , TRIM(myName)//'_t'//clstr                 &
                           , p4=Par(ic)%CPLDATA(ii)%target                   &
                           , reprid = reprid                                 &
                           , lrestreq=.TRUE.)
                      CALL channel_halt (substr, status)
                      Par(ic)%CPLData(ii)%ltimedep = .FALSE.

                      IF (Par(ic)%CPLDATA(ii)%l_sendunit) THEN
                         CALL MMD_Inter_Bcast(unit,.FALSE.,child_id)

                         CALL new_attribute(status, submodstr &
                              , TRIM(myName)//'_t'//clstr , 'units' &
                              , c=unit)
                         CALL channel_halt(substr, status)
                      END IF
                   ENDIF if_mmd2way
#ifdef ECHAM5
                END IF if_tracer
#endif
                ! Make memory for data exchange field
                ! 1. get representation id. Use namelist entry
                !CALL get_representation_id(status &
                !     , Par(ic)%CPLData(ii)%repr, reprid)
                !CALL channel_halt (substr, status)
                ! 3. make the memory
                IF (Par(ic)%CPLData(ii)%repr == 'GP_2D_HORIZONTAL') THEN
                   reprid = GP_2D_HORIZONTAL
                ELSE IF (Par(ic)%CPLData(ii)%repr == 'GP_3D_MID') THEN
                   reprid = GP_3D_RDCX(ic)
                ELSE IF (Par(ic)%CPLData(ii)%repr == 'GP_3D_INT') THEN
                   reprid = GP_3D_RDCXp1(ic)
                END IF

                CALL int2str(clstr, child_id)
                CALL new_channel_object(status, submodstr           &
                     , TRIM(myName)//'_'//clstr                     &
                     , p4=Par(ic)%CPLDATA(ii)%ptr, reprid = reprid  &
                     , lrestreq=.TRUE.)
                CALL channel_halt (substr, status)
                ! 4. get additional information: axis string and
                !    local dimensions
                CALL get_representation_info(status, ' ', reprid &
                     , rank=Par(ic)%CPLDATA(ii)%rank             &
                     , axis=Par(ic)%CPLDATA(ii)%AXIS             &
                     , gdimlen=Par(ic)%CPLDATA(ii)%gdimlen       &
                     , ldimlen=Par(ic)%CPLDATA(ii)%ldimlen)
                CALL channel_halt(substr, status)

             ! Define exchange field and its dimensions in library
             CALL MMD_P_Set_ParDataArray(ic, status    &
                  , Par(ic)%CPLDATA(ii)%ldimlen        &
                  , Par(ic)%CPLDATA(ii)%axis           &
                  , p4=Par(ic)%CPLDATA(ii)%ptr)
             IF (status /=0)  &
                  CALL error_bi('ERROR in MMD_P_Set_ParDataArray',substr)

#ifdef ECHAM5
             !             IF (Par(ic)%lgrhum) THEN
             ! vapour?
             IF (TRIM(myName) == 'qm1')  Par(ic)%idx_qv = ii
             ! cloud water?
             IF (TRIM(myName) == 'xlm1') Par(ic)%idx_qc = ii
             ! cloud ice?
             IF (TRIM(myName) == 'xim1') Par(ic)%idx_qi = ii
             !             ENDIF
             ! horiz. wind component ?
             IF (TRIM(myName) == 'u' .OR. TRIM(myName) == 'um1')       &
                  Par(ic)%idx_u    = ii
             IF (TRIM(myName) == 'v' .OR. TRIM(myName) == 'vm1')       &
                  Par(ic)%idx_v    = ii
             IF (TRIM(myName) == 't' .OR. TRIM(myName) == 'tm1')       &
                  Par(ic)%idx_t    = ii
             IF (TRIM(myName) == 'alps' .OR. TRIM(myName) == 'alpsm1') &
                  Par(ic)%idx_ps   = ii
#endif
#ifdef COSMO
             IF (Par(ic)%lgrhum) THEN
                ! vapour?
                IF (TRIM(myName) == 'QV') Par(ic)%idx_qv = ii
                ! cloud water?
                IF (TRIM(myName) == 'QC') Par(ic)%idx_qc = ii
                ! cloud ice?
                IF (TRIM(myName) == 'QI') Par(ic)%idx_qi = ii
             END IF
#endif
          ENDDO ! END WHILE
       ENDIF ! NEXCH

       ! test if required fields for calulation of generalised humidity
       ! are available
       IF (Par(ic)%lgrhum) THEN
          IF (Par(ic)%idx_qv == -99 .OR. Par(ic)%idx_qc == -99) THEN
             CALL error_bi(&
             'FIELDs required for generalised humidity coupling not available'&
                  , substr)
          END IF
       ENDIF
    END DO ! END CHILD LOOP

    CALL end_message_bi(modstr, 'INIT Parent COUPLING', substr)
    ! PARENT COUPLING -
    ! ---------------------------------------------------------------------
    ENDIF! ldev
    CALL start_message_bi(modstr, 'ALLOCATE BUFFER(S)', substr)

    DO ic=1,NumberOfChildren
       CALL MMD_P_SetInd_and_AllocMem(ic)
    ENDDO

    CALL end_message_bi(modstr, 'ALLOCATE BUFFER(S)', substr)

    RETURN

  END SUBROUTINE mmd2way_parent_init_coupling
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE mmd2way_parent_global_start

    ! MESSy/SMCL
    USE messy_main_timer,         ONLY:  lfirst_cycle
    ! MMD
    USE mmd_parent,               ONLY: MMD_P_FillBuffer
    IMPLICIT NONE

    ! I/O

    ! LOCAL
    INTEGER                     :: ic
    CHARACTER(LEN=*), PARAMETER :: substr='mmd2way_parent_global_start'

    do_childs: DO ic=1,NumberOfChildren

       CL(ic)%lcpl = CL(ic)%lcpl .OR. lfirst_cycle
       IF (.NOT.ldev) CL(ic)%waittime = 0.0_dp ! ub_ak_20170324

       ifcpl: IF (CL(ic)%lcpl ) THEN

#ifdef COSMO
          CALL debug_bi(modstr, 'COSMO PARENT START EXCHANGE DATA','',substr)
#endif
#ifdef ECHAM5
          CALL debug_bi(modstr, 'ECHAM PARENT START EXCHANGE DATA','',substr)
#endif
          ! Send Initial Data to child
          IF (ldev) THEN
             CALL MMD_P_FillBuffer (ic)
          ELSE
             CALL MMD_P_FillBuffer (ic, WaitTime = Cl(ic)%Waittime)
          END IF
#ifdef COSMO
          CALL debug_bi(modstr, 'COSMO PARENT STOP EXCHANGE DATA','', substr)
#endif
#ifdef ECHAM5
          CALL debug_bi(modstr, 'ECHAM PARENT STOP EXCHANGE DATA','', substr)
#endif
      ENDIF ifcpl

      ! ---------------------------------------------------------------------
    END DO do_childs

    RETURN

  END SUBROUTINE mmd2way_parent_global_start
! ---------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE mmd2way_parent_global_end(flag)

    ! MESSy/BMIL
    USE messy_main_timer_bi,      ONLY: event_state
    ! MESSy/SMCL
    USE messy_main_timer,         ONLY: lbreak, l_rerun, lstop &
                                      , next_date
    ! MMD
    USE mmd_parent,               ONLY: MMD_P_GetBuffer      &
                                      , MMD_Parent_for_Child &
                                      , MMD_Inter_Bcast

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)         :: flag !only important for COSMO-COSMO coupling

    ! LOCAL
    INTEGER                     :: ic
    INTEGER,DIMENSION(3)        :: Timeflags
    INTEGER                     :: child_id
    CHARACTER(LEN=*), PARAMETER :: substr='mmd2way_parent_global_end'

    do_childs: DO ic=1,NumberOfChildren
       child_id = MMD_Parent_for_Child(ic)

       SELECT CASE(flag)

       CASE (1,3)

          IF (flag ==1) THEN
             CL(ic)%lcpl = event_state(CPL_EVENT(ic), next_date)
             IF (ldev) CL(ic)%waittime = 0._dp
          END IF
          IF (.NOT. ldev) THEN
             RETURN
          END IF
          ifcpl: IF (CL(ic)%lcpl ) THEN

          ! -------------------------------------------------------------------
          ! PARENT coupling+

          CALL debug_bi(modstr, 'PARENT START EXCHANGE DATA','',substr)

          IF (flag == 1) THEN

             CALL MMD_P_GetBuffer (ic, WaitTime = CL(ic)%Waittime)

          ENDIF

          CALL mmd2way_parent_couple_gp(ic, flag)

          CALL debug_bi(modstr, 'PARENT STOP EXCHANGE DATA','',substr)
          ! PARENT coupling-
          ! -------------------------------------------------------------------

       ENDIF ifcpl

       ! ---------------------------------------------------------------------
    CASE (2)
          CALL debug_bi(modstr, 'START EXCHANGE BREAK INFO','', substr)

          ! Send job control information to child
          Timeflags = 0
          IF (lbreak)    Timeflags(1) = 1
          IF (l_rerun)   Timeflags(2) = 1
          IF (lstop)     Timeflags(3) = 1

          CALL MMD_Inter_Bcast (Timeflags, .TRUE., child_id)
          CALL debug_bi(modstr, 'END EXCHANGE BREAK INFO','', substr)

       ! ---------------------------------------------------------------------
       END SELECT
    END DO do_childs

    RETURN

  END SUBROUTINE mmd2way_parent_global_end
  ! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE mmd2way_parent_free_memory

    USE mmd_parent,  ONLY: MMD_P_FreeMem, MMD_Parent_for_Child
    USE mmd_test,    ONLY: MMD_testP_FreeMem

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    INTEGER :: ic, child_id

    DO ic =1, NumberofChildren
       IF (ldev .AND. ALLOCATED(PAR)) THEN
          IF (ASSOCIATED(PAR(ic)%all_plon)) THEN
             DEALLOCATE(PAR(ic)%all_plon)
             PAR(ic)%all_plon => NULL()
          ENDIF
          IF (ASSOCIATED(PAR(ic)%all_plat)) THEN
             DEALLOCATE(PAR(ic)%all_plat)
             PAR(ic)%all_plat => NULL()
          ENDIF
          IF (ASSOCIATED(PAR(ic)%CPLDATA)) DEALLOCATE(PAR(ic)%CPLDATA)
          DEALLOCATE(PAR)
       ENDIF
       child_id = MMD_Parent_for_Child(ic)
       CALL MMD_testP_FreeMem(child_id)
       IF (ASSOCIATED(CL(ic)%CPLDATA)) DEALLOCATE(CL(ic)%CPLDATA)
    END DO

    CALL MMD_P_FreeMem(NumberofChildren)

    DEALLOCATE(CL)
    DEALLOCATE(CPL_IOEVENT)
    DEALLOCATE(CPL_EVENT)

  END SUBROUTINE mmd2way_parent_free_memory
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! Private Subroutines
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

  SUBROUTINE Setup_Child_Timer (child_id, id)

    ! MESSy/BMIL
    USE messy_main_mpi_bi,        ONLY: p_pe, p_io, p_bcast
    USE messy_main_timer_bi,      ONLY: timer_event_init
    ! MESSy/SMCL
    USE messy_main_timer,         ONLY: current_date, start_date &
                                      , stop_date,   resume_date &
                                      , delta_time

    ! MMD
    USE mmd_parent,               ONLY: MMD_Recv_from_Child &
                                      , MMD_Inter_Bcast

    IMPLICIT NONE

    INTRINSIC :: INT

    INTEGER,INTENT(IN)            :: child_id ! Child id (all childs)
    INTEGER,INTENT(IN)            :: id       ! Child ID (this parent)

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER   :: substr ='mmd2way_parent_SetupChildTimer'
    INTEGER                       :: status
    INTEGER,DIMENSION(9)          :: timing_array1=0
    INTEGER,DIMENSION(1)          :: cplseconds


    IF (p_pe == 0)  THEN
       CALL MMD_Recv_from_Child (child_id, cplseconds,1, 0, 61, status)
       IF (status /= 0) CALL error_bi('MMD_Recv_from_Child ERROR', substr)
    ENDIF
    CALL p_bcast(cplseconds, p_io)

    ! INITIAL COUPLER TIMER FOR EACH CHILD
    CPL_IOEVENT(id) = io_time_event(CPLSECONDS(1), 'seconds', 'first', 0)

    ! coupling evaluation moved to global end
    CALL timer_event_init( CPL_EVENT(id), CPL_IOEVENT(id) &
                         , 'CPLPAR', 'next')

    ! Send dates to child for synchonization
    ! - current_date
    timing_array1(1) = current_date%day
    timing_array1(2) = current_date%second

    ! - start_date
    timing_array1(3) = start_date%day
    timing_array1(4) = start_date%second
    ! - resume_date
    timing_array1(5) = resume_date%day
    timing_array1(6) = resume_date%second
    ! - stop_date
    timing_array1(7) = stop_date%day
    timing_array1(8) = stop_date%second
    ! - time step
    timing_array1(9) = INT(delta_time)

    CALL MMD_Inter_Bcast (timing_array1, .TRUE., child_id)

    RETURN

  END SUBROUTINE Setup_Child_Timer
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE Setup_Child_Area (child_id, ic)

    ! MESSy/BMIL
    USE messy_main_mpi_bi,        ONLY: p_pe, p_io, p_bcast
#ifdef ECHAM5
    USE messy_main_transform_bi,  ONLY: locate_in_decomp
    USE messy_main_grid_def_mem_bi, ONLY: vct, ngl, nlev
    USE messy_main_grid_def_bi,     ONLY: philat, philon
#endif
#ifdef COSMO
    USE messy_main_grid_def_mem_bi, ONLY: startlon_tot, startlat_tot   &
                                      , pollon, pollat, polgam       &
                                      , dlon, dlat, ke_tot           &
                                      , ie_tot, je_tot               &
                                      , vcoord, refatm               &
                                      , svc1, svc2                   &
                                      , ke_soil
    USE messy_main_grid_def_bi,     ONLY: czmls
#endif

    ! MMD
    USE mmd_parent,               ONLY: MMD_Send_to_Child   &
                                      , MMD_Recv_from_Child
#ifdef ECHAM5
    ! MESSy/SMCL
    USE messy_main_constants_mem, ONLY: pi
#endif

    IMPLICIT NONE

    INTRINSIC :: REAL, SIZE
#ifdef ECHAM5
    INTRINSIC :: ABS, MAX, MAXVAL, MIN, MINVAL
#endif
    INTEGER,INTENT(IN)  :: child_id ! Child id (all childs)
    INTEGER,INTENT(IN)  :: ic       ! Child id (parent specific)

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'Setup_Child_Area'

    INTEGER, DIMENSION(3)                      :: val
    INTEGER                                    :: status
    INTEGER                                    :: ie_tot_cl,je_tot_cl
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: rlon_glob,rlat_glob
    REAL(KIND=DP), DIMENSION(17)               :: Define_coarse_grid_real
    INTEGER, DIMENSION(8)                      :: Define_coarse_grid_int
#ifdef COSMO
    REAL(KIND=DP), PARAMETER                   :: escal = 3._dp
#endif

#ifdef ECHAM5
    INTEGER, DIMENSION(4,2)                    :: corner
    INTEGER                                    :: i, j
    REAL(KIND=DP)                              :: zminlon
    REAL(KIND=DP)                              :: zmaxlon
    REAL(KIND=DP)                              :: zminlat
    REAL(KIND=DP)                              :: zmaxlat
    REAL(KIND=DP)                              :: startlon_exch
    INTEGER                                    :: zpe, zjp, zjrow
    REAL(KIND=DP)                              :: dlon, dlat
    LOGICAL                                    :: l_pacific
    LOGICAL                                    :: l_mixedsign
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: rlon_glob_pos,rlon_glob_neg

    REAL(KIND=DP), PARAMETER, DIMENSION(5) :: soil_depth &
         = (/0.03_dp, 0.19_dp, 0.78_dp, 2.68_dp, 6.98_dp/)
    ! NBOUNDLINES(INT2COSMO) + 2.5  = 3.5 ???
    ! REAL(KIND=DP), PARAMETER :: escal = 4.5
    ! nzbounds in INT2LM is set to 2 and nboundlines is 1,
    ! additionally we need some extra space for adjustment within the
    ! box. Therefore 4.0 should be a secure estimate for the frame of the
    ! boundary data
    REAL(KIND=DP), PARAMETER               :: escal = 4._dp
#endif

    !
    ! =======================================================================

    IF (p_pe == 0) THEN
      CALL MMD_Recv_from_Child (child_id, val, SIZE(val), 0, 123, status)
      IF (status /= 0) CALL error_bi('MMD_Recv_from_Child ERROR 01', substr)

      ie_tot_cl   = val(1)
      je_tot_cl   = val(2)
      CL(ic)%NEXCH = val(3)

      ! Get geographical coord from child
      ALLOCATE(rlat_glob(ie_tot_cl,je_tot_cl))
      ALLOCATE(rlon_glob(ie_tot_cl,je_tot_cl))

      CALL MMD_Recv_from_Child ( child_id, rlat_glob &
                                , SIZE(rlat_glob), 0, 11, status)
      IF (status /= 0) CALL error_bi('MMD_Recv_from_Child ERROR 02', substr)
      CALL MMD_Recv_from_Child ( child_id, rlon_glob  &
                                , SIZE(rlon_glob), 0, 12, status)
      IF (status /= 0) CALL error_bi('MMD_Recv_from_Child ERROR 03', substr)

#ifdef COSMO
      ! IMPLEMENTATION FOR COSMO PARENT

      ! Set up Coarse grid information for child.
      ! 1. version, send complete Parent area

      Define_coarse_grid_real(1) = startlat_tot    ! Latitude  South-West
      Define_coarse_grid_real(2) = startlon_tot    ! Longitude South-West
      Define_coarse_grid_real(3) = pollat
      Define_coarse_grid_real(4) = pollon
      Define_coarse_grid_real(5) = polgam
      Define_coarse_grid_real(6) = dlat
      Define_coarse_grid_real(7) = dlon
      ! Latitude  North-East
      Define_coarse_grid_real(8) = startlat_tot+je_tot*dlat
      ! Longitude North-East
      Define_coarse_grid_real(9) = startlon_tot+ie_tot*dlon
      ! vertical grid
      Define_coarse_grid_real(10) = vcoord%vcflat
      Define_coarse_grid_real(11) = refatm%p0sl
      Define_coarse_grid_real(12) = refatm%t0sl
      Define_coarse_grid_real(13) = refatm%dt0lp
      Define_coarse_grid_real(14) = refatm%delta_t
      Define_coarse_grid_real(15) = refatm%h_scal
      Define_coarse_grid_real(16) = svc1
      Define_coarse_grid_real(17) = svc2

      Define_coarse_grid_int(1)  = ie_tot
      Define_coarse_grid_int(2)  = je_tot
      Define_coarse_grid_int(3)  = ke_tot
      Define_coarse_grid_int(4)  = vcoord%ivctype
      Define_coarse_grid_int(5)  = refatm%irefatm

      Define_coarse_grid_int(6)  = ke_soil
      Define_coarse_grid_int(7)  = 2      ! dummy for itype_w_so_rel (needed
                                          ! if ECHAM is parent)
      Define_coarse_grid_int(8)  = 0      ! dummy for itype_t_cl     (needed
                                          ! if ECHAM is parent)
#endif

#ifdef ECHAM5
! IMPLEMENTATION FOR ECHAM PARENT

      ! CHECK IF COSMO DOMAIN CONTAINS 0 deg LONGITUDE
      l_pacific = ( (rlon_glob(1,1) > rlon_glob(ie_tot_cl,1)) .OR. &
                    (rlon_glob(1,je_tot_cl) > rlon_glob(ie_tot_cl,je_tot_cl)))

      IF (l_pacific) THEN
         ALLOCATE(rlon_glob_pos(ie_tot_cl,je_tot_cl))
         ALLOCATE(rlon_glob_neg(ie_tot_cl,je_tot_cl))
         DO i = 1,ie_tot_cl
            DO j = 1,je_tot_cl
               IF (rlon_glob(i,j) > 0._dp) THEN
                  rlon_glob_pos(i,j) =  rlon_glob(i,j)
                  rlon_glob_neg(i,j) = -9999.0_dp
               ELSE
                  rlon_glob_neg(i,j) =  rlon_glob(i,j)
                  rlon_glob_pos(i,j) = 9999.0_dp
               ENDIF
            ENDDO
         ENDDO
         zminlon = MINVAL(rlon_glob_pos(:,:))*180.0_dp/pi
         zminlat = MINVAL(rlat_glob(:,:))*180.0_dp/pi
         zmaxlon = MAXVAL(rlon_glob_neg(:,:))*180.0_dp/pi
         zmaxlat = MAXVAL(rlat_glob(:,:))*180.0_dp/pi
         DEALLOCATE (rlon_glob_pos)
         DEALLOCATE (rlon_glob_neg)
      ELSE
         zminlon = MINVAL(rlon_glob(:,:))*180.0_dp/pi
         zminlat = MINVAL(rlat_glob(:,:))*180.0_dp/pi
         zmaxlon = MAXVAL(rlon_glob(:,:))*180.0_dp/pi
         zmaxlat = MAXVAL(rlat_glob(:,:))*180.0_dp/pi
      ENDIF

      ! ADD BOUNDARY FOR INTERPOLATION
      dlon = 360.0_dp/REAL(SIZE(philon),dp)
      dlat = (MAXVAL(philat) - MINVAL(philat))/REAL(ngl-1, dp)
      zminlon = zminlon - escal* dlon
      IF (zminlon < -180.0_dp) zminlon = zminlon + 360.0_dp
      zminlat = MAX(zminlat - escal * dlat, MINVAL(philat))
      zmaxlon = zmaxlon + escal * dlon
      IF (zmaxlon >  180.0_dp) zmaxlon = zmaxlon - 360.0_dp
      zmaxlat = MIN(zmaxlat + escal* dlat, MAXVAL(philat))

      l_mixedsign = zminlon < 0._dp .AND. zmaxlon > 0._dp

      !write(0,*) 'INT2COSMO DOMAIN:', zminlon,zmaxlon,zminlat,zmaxlat

      CALL locate_in_decomp(status, zminlon, zminlat, zpe, zjp, zjrow &
           , jgx = corner(1,1), jgy = corner(1,2) )
      IF (status /= 0) CALL error_bi(&
           'ERROR in locate_in_decomp (1)','mmd2way_parent')

      CALL locate_in_decomp(status, zminlon, zmaxlat, zpe, zjp, zjrow &
           , jgx = corner(2,1), jgy = corner(2,2) )
      IF (status /= 0) CALL error_bi(&
           'ERROR in locate_in_decomp (2)','mmd2way_parent')

      CALL locate_in_decomp(status, zmaxlon, zmaxlat, zpe, zjp, zjrow &
           , jgx = corner(3,1), jgy = corner(3,2) )
      IF (status /= 0) CALL error_bi(&
           'ERROR in locate_in_decomp (3)','mmd2way_parent')

      CALL locate_in_decomp(status, zmaxlon, zminlat, zpe, zjp, zjrow &
           , jgx = corner(4,1), jgy = corner(4,2) )
      IF (status /= 0) CALL error_bi(&
           'ERROR in locate_in_decomp (4)','mmd2way_parent')

      ! recalculate dlat (gaussian grid is not equidistant!)
      dlat = ABS((philat(corner(2,2)) - philat(corner(1,2))) /  &
           REAL(corner(2,2)-corner(1,2), dp) )

      ! Set up Coarse grid information for child.
      ! Add an additional box in every direction to the ECHAM5 space
      IF (philon(corner(1,1)) > 180._dp) THEN
         startlon_exch = philon(corner(1,1)) - 360._dp
      ELSE
         startlon_exch = philon(corner(1,1))
      ENDIF
      Define_coarse_grid_real(1) = philat(corner(1,2))    ! startlat_tot
      Define_coarse_grid_real(2) = startlon_exch          ! startlon_tot
      Define_coarse_grid_real(3) = 90.0_dp                ! pollat
      Define_coarse_grid_real(4) = 180.0_dp               ! pollon
      Define_coarse_grid_real(5) = 0.0                    ! polgam
      Define_coarse_grid_real(6) = dlat                   ! dlat
      Define_coarse_grid_real(7) = dlon                   ! dlon
      Define_coarse_grid_real(8) = philat(corner(3,2))    ! endlat_tot
      Define_coarse_grid_real(9) = philon(corner(3,1))    ! endlon_tot
      ! DUMMIES (used for vertical grid COSMO llm2lm)
      Define_coarse_grid_real(10) = 0.
      Define_coarse_grid_real(11) = 0.
      Define_coarse_grid_real(12) = 0.
      Define_coarse_grid_real(13) = 0.
      Define_coarse_grid_real(14) = 0.
      Define_coarse_grid_real(15) = 0.
      Define_coarse_grid_real(16) = 0.
      Define_coarse_grid_real(17) = 0.

      IF ( l_mixedsign .AND. (.NOT.l_pacific) ) THEN
         Define_coarse_grid_int(1)  = &                   ! ie_tot
              SIZE(philon) - ( corner(1,1) - corner(4,1) + 1) + 2
      ELSE
         Define_coarse_grid_int(1)  = &                   ! ie_tot
              corner(4,1) - corner(1,1) + 1
      END IF

      Define_coarse_grid_int(2)  = (corner(1,2))-(corner(2,2))+1 ! je_tot
      Define_coarse_grid_int(3)  = nlev                          ! ke_tot
      Define_coarse_grid_int(4)  = 0        ! dummy for ivctype in llm2lm
      Define_coarse_grid_int(5)  = 0        ! dummy for irefatm in llm2lm
      Define_coarse_grid_int(6)  = 4        ! number of soil layers
      Define_coarse_grid_int(7)  = 2        ! itype_w_so_rel
      Define_coarse_grid_int(8)  = 1        ! itype_t_cl

#endif
      ! Send coarse grid information to child
      CALL MMD_Send_to_Child (child_id, Define_coarse_grid_real &
                              , 17, 0, 21, status)
      IF (status /= 0) CALL error_bi('MMD_Send_to_Child ERROR 01', substr)

      CALL MMD_Send_to_Child (child_id, Define_coarse_grid_int  &
                              ,  8, 0, 22, status)
      IF (status /= 0) CALL error_bi('MMD_Send_to_Child ERROR 02', substr)
#ifdef ECHAM5
      CALL MMD_Send_to_Child (child_id, vct,  2*(nlev+1), 0, 23, status)

      IF (status /= 0) CALL error_bi('MMD_Send_to_Child ERROR 03 E', substr)
      CALL MMD_Send_to_Child (child_id, soil_depth,  5, 0, 26, status)
      IF (status /= 0) CALL error_bi('MMD_Send_to_Child ERROR 04 E', substr)
#endif
#ifdef COSMO
      ! Send vcoord to child (analog to vct in ECHAM)
      IF (vcoord%ivctype == 1) THEN
         CALL MMD_Send_to_Child (child_id, vcoord%sigm_coord &
              , INT(ke_tot+1), 0, 23, status)
      ELSE
         CALL MMD_Send_to_Child (child_id, vcoord%vert_coord &
              , INT(ke_tot+1), 0, 23, status)
      END IF

      IF (status /= 0) CALL error_bi('MMD_Send_to_Child ERROR 03 C', substr)
      ! Send depth of soil layers
      CALL MMD_Send_to_Child (child_id, czmls, INT(ke_soil+1), 0, 26, status)
      IF (status /= 0) CALL error_bi('MMD_Send_to_Child ERROR 04 C', substr)
#endif

      CALL send_lat_and_lon

      DEALLOCATE (rlon_glob)
      DEALLOCATE (rlat_glob)
   ELSE
#ifdef COSMO
      !     No dummy CALL on all PE's to do initialization phase
#endif
#ifdef ECHAM5
      ! One dummy call on all PE's to do initialization phase
      CALL locate_in_decomp (status, 0.0_dp, 0.0_dp, zpe, zjp, zjrow &
           ,jgx = corner(1,1), jgy = corner(1,2))
#endif
   END IF

   ! BROADCAST NUMBER OF COUPLING ARRAYs
   CALL p_bcast(CL(ic)%NEXCH, p_io)

   RETURN

 CONTAINS

!--------------------------------------------------------------------------
#ifdef COSMO
    SUBROUTINE send_lat_and_lon

      IMPLICIT NONE

      INTRINSIC :: NINT

      ! LOCAL
      INTEGER                      :: i,j,ian,jan
      INTEGER                      :: idlat, idlon, istartlon, istartlat
      REAL(dp), DIMENSION(je_tot)  :: lats
      REAL(dp), DIMENSION(ie_tot)  :: lons
      REAL(dp), PARAMETER          :: fac = 100000._dp

      jan = je_tot
      ian = ie_tot

      ! bugfix replace INT by NINT to avoid rounding errors
      istartlon = NINT(startlon_tot * fac)
      istartlat = NINT(startlat_tot * fac)
      idlon = NINT(dlon * fac)
      idlat = NINT(dlat * fac)

      DO j=1,jan
        lats(j) = REAL(istartlat+(j-1)*idlat, dp)/ fac
      END DO

      DO i=1,ian
        lons(i) = REAL(istartlon + (i-1)*idlon, dp)/ fac
      END DO

      ! ... and send to child
      CALL MMD_Send_to_Child (child_id, lats, jan, 0, 24, status)
      IF (status /= 0) CALL error_bi('MMD_Send_to_Child ERROR 01', 'lat / lon')
      CALL MMD_Send_to_Child (child_id, lons, ian, 0, 25, status)
      IF (status /= 0) CALL error_bi('MMD_Send_to_Child ERROR 02', 'lat / lon')

      RETURN

   END SUBROUTINE send_lat_and_lon
!--------------------------------------------------------------------------
#endif
#ifdef ECHAM5
!--------------------------------------------------------------------------

    SUBROUTINE send_lat_and_lon

      IMPLICIT NONE

      ! LOCAL
      INTEGER                                      :: i,j,ian,jan
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE     :: lats
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE     :: lons

      jan =  Define_coarse_grid_int(2)
      ALLOCATE(lats(jan))

      ian = Define_coarse_grid_int(1)
      ALLOCATE(lons(ian))

      ! COPY lats
      DO j=1,jan
         lats(j) = philat(corner(1,2) - j + 1)
      END DO

      ! FLIP lons    [0 ... 360] -> [-180 ... 180]
      DO i=1,ian
         IF (corner(1,1) + i - 1 <= SIZE(philon)) THEN
            lons(i) = philon(corner(1,1) + i - 1)
         ELSE
            lons(i) = philon(corner(1,1) + i - 1 - SIZE(philon))
         END IF
         IF (lons(i) > 180.0_dp) lons(i) = lons(i) - 360.0_dp
      END DO

      !... and send to child
      CALL MMD_Send_to_Child (child_id, lats, jan, 0, 24, status)
      IF (status /= 0) &
           CALL error_bi('MMD_Send_to_Child ERROR 01 E', 'lat / lon')
      CALL MMD_Send_to_Child (child_id, lons, ian, 0, 25, status)
      IF (status /= 0) &
           CALL error_bi('MMD_Send_to_Child ERROR 02 E', 'lat / lon')

      DEALLOCATE(lats)
      DEALLOCATE(lons)

      RETURN

   END SUBROUTINE send_lat_and_lon
!--------------------------------------------------------------------------
#endif

#ifdef CESM1
   SUBROUTINE send_lat_and_lon
   END SUBROUTINE send_lat_and_lon
#endif

  END SUBROUTINE Setup_Child_Area
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
  SUBROUTINE Setup_Data_Exchange_with_Child (child_id,ic)

    ! MESSy/BMIL
    USE messy_main_transform_bi,  ONLY: locate_in_decomp
    USE messy_main_mpi_bi,        ONLY: p_pe
#ifdef ECHAM5
    USE messy_main_mpi_bi,        ONLY: dcl
#endif
    ! MMD
    USE mmd_parent,               ONLY: MMD_Recv_from_Child  &
                                      , MMD_P_Set_Indexlist
    USE mmd_test,                 ONLY: MMD_testP_Setup      &
                                      , MMD_testP_Fill       &
                                      , MMD_testP_FinishFill
#ifdef COSMO
    ! COSMO
    USE data_modelconfig,         ONLY: ie_max,je_max &
                                     , startlon_tot, startlat_tot   &
                                     , startlon,     startlat       &
                                     , endlon_tot, endlat_tot
#endif

    IMPLICIT NONE

    INTRINSIC :: SIZE

    INTEGER,INTENT(IN)    :: ic
    INTEGER,INTENT(IN)    :: child_id

    ! LOCAL
    INTEGER               :: i, j, k, status
    INTEGER               :: icount                         !Index cout
    INTEGER, DIMENSION(3) :: size_of_array
    INTEGER, DIMENSION(5) :: cgrid

    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: all_lon,all_lat
    INTEGER, DIMENSION(:,:), ALLOCATABLE         :: index_list

#ifdef COSMO
    CHARACTER(LEN=*), PARAMETER :: substr='cosmo_parent_exchange'
    CALL MMD_testP_Setup (child_id,INT(ie_max),INT(je_max),'XYNZ')
#endif
#ifdef ECHAM5
    CHARACTER(LEN=*), PARAMETER :: substr='echam_parent_exchange'
    CALL MMD_testP_Setup (child_id, dcl%nproma, dcl%ngpblks,'XZNY')
#endif
    IF (p_pe == 0) THEN
      CALL MMD_Recv_from_Child (child_id, size_of_array, 3, 0, 40, status)
      IF (status /= 0) CALL error_bi('MMD_Recv_from_Child 01',substr)

      ALLOCATE (all_lon(size_of_array(1),size_of_array(2),size_of_array(3)))
      ALLOCATE (all_lat(size_of_array(1),size_of_array(2),size_of_array(3)))

      CALL MMD_Recv_from_Child(child_id, all_lon, SIZE(all_lon), 0, 41, status)
      IF (status /= 0) CALL error_bi('MMD_Recv_from_Child 02',substr)
      CALL MMD_Recv_from_Child(child_id, all_lat, SIZE(all_lat), 0, 42, status)
      IF (status /= 0) CALL error_bi('MMD_Recv_from_Child 03',substr)

      ALLOCATE (index_list(6, SIZE(all_lon)))

      ! Create list of required PARENT grid boxes
      icount = 0
      DO k=1,size_of_array(3)
         DO j=1,size_of_array(2)
            DO i=1,size_of_array(1)
               ! Legal Value
               if_legal: IF (all_lon(i,j,k) > -1000._dp .AND. &
                    all_lat(i,j,k) > -1000._dp) THEN

                  CALL locate_in_decomp (status, all_lon(i,j,k), all_lat(i,j,k) &
                      ,cgrid(1),cgrid(2),cgrid(3)                               &
#ifdef COSMO
                      ,lrot = .TRUE.                                            &
#endif
                      )
                  IF (status /= 0) THEN
                     write(iouerr,*) ' locate_in_decomp', i,j,k &
                          ,all_lon(i,j,k), all_lat(i,j,k), 'LON MIN/MAX k ',k &
                          , MINVAL( all_lon(:,:,k)), MAXVAL( all_lon(:,:,k))  &
                          , 'LON MIN/MAX ALL ',MINVAL( all_lon(:,:,:))        &
                          , MAXVAL( all_lon(:,:,:)), 'LAT MIN/MAX k ',k       &
                          , MINVAL( all_lat(:,:,k)), MAXVAL( all_lat(:,:,k))  &
                          , 'LAT MIN/MAX ALL ',MINVAL( all_lat(:,:,:))        &
                          , MAXVAL( all_lat(:,:,:))
#ifdef COSMO
                     write(iouerr,*) ' locate_in_decomp', startlon_tot, endlon_tot &
                          , startlon, 'LAT ',  startlat_tot, endlat_tot, startlat
#endif
                     CALL error_bi('locate in decomp',substr)
                  ENDIF
                  !write(0,'(a,3i5,a,5i5,a,i5,2f7.2)') 'Locate ',i,j,icount+1  &
                  !     ,' cgrid ',cgrid, ' COSMO_PE ',k-1, all_lon(i,j,k) &
                  !     , all_lat(i,j,k)
                  icount = icount+1
                  index_list(1,icount) = cgrid(2) ! First index in Parent Array
                  index_list(2,icount) = cgrid(3) ! Second index in Parent Array
                  index_list(3,icount) = i        ! Lon Index child coarse grid
                  index_list(4,icount) = j        ! Lat Index child coarse grid
                  index_list(5,icount) = k-1      ! PE Number child
                  index_list(6,icount) = cgrid(1) ! PE Number parent
                  CALL MMD_testP_Fill (child_id, index_list(:,icount) &
                       ,  all_lon(i,j,k), all_lat(i,j,k))
               END IF if_legal
            END DO
         END DO
      END DO
      CALL MMD_P_Set_Indexlist (ic,index_list(:,1:icount))

      DEALLOCATE(all_lon)
      DEALLOCATE(all_lat)
    ELSE
      ALLOCATE (index_list(6, 1))                             ! dummy allocate
      CALL MMD_P_Set_Indexlist (ic,index_list)
    END IF

    CALL MMD_testP_FinishFill (child_id)

    DEALLOCATE(index_list)

    RETURN

  END SUBROUTINE Setup_Data_Exchange_with_Child
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE Define_data_arrays (child_id, ic)

    ! MESSy/BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! MESSy/SMCL
    USE messy_main_channel,          ONLY: STRLEN_CHANNEL,STRLEN_OBJECT &
                                         , get_channel_object           &
                                         , get_channel_object_info      &
                                         , get_attribute                &
                                         , set_channel_object_restreq   &
                                         , get_channel_object_dimvalue
    USE messy_main_tools,            ONLY: str
    USE messy_main_channel_repr,     ONLY: get_representation_info
    USE messy_main_constants_mem,    ONLY: STRLEN_MEDIUM, STRLEN_ULONG
    ! MMD
    USE mmd_parent,               ONLY: MMD_P_GetNextArray           &
                                      , MMD_P_Set_DataArray          &
                                      , MMD_P_Send_repr
    USE mmd_test,                 ONLY: MMD_testP_GetTestPtr
    USE mmd_utilities,            ONLY: STRLEN_ATT
    USE mmd_mpi_wrapper,          ONLY: MMD_Inter_Bcast

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, TRIM

    INTEGER, INTENT(IN)           :: child_id
    INTEGER, INTENT(IN)           :: ic
    ! LOCAL
    INTEGER                       :: ii, ix
    INTEGER                       :: status
    CHARACTER(LEN=STRLEN_CHANNEL) :: myChannel
    CHARACTER(LEN=STRLEN_OBJECT)  :: myName
    CHARACTER(LEN=STRLEN_MEDIUM)  :: clntrepr
    INTEGER                       :: reprid
    CHARACTER(LEN=*), PARAMETER   :: substr ='mmd2way_parent:Define_data_arrays'
    CHARACTER(LEN=STRLEN_MEDIUM)  :: repr_name
    INTEGER, DIMENSION(4)         :: repr_gdimlen
    CHARACTER(LEN=STRLEN_ATT)     :: ch_att = ' '
    LOGICAL                       :: l_sendunit
    CHARACTER(LEN=STRLEN_ULONG)   :: unit    = ' '
    CHARACTER(LEN=STRLEN_ULONG)   :: levtype = ' '
    REAL(dp), DIMENSION(:), POINTER :: levs => NULL()

    ALLOCATE (CL(ic)%CPLDATA(CL(ic)%NEXCH))

    ii=0
    DO WHILE (MMD_P_GetNextArray (ic, myChannel, myName, clntrepr, l_sendunit))
       ii=ii+1

       CL(ic)%CPLDATA(ii)%NAME%CHA=myChannel
       CL(ic)%CPLDATA(ii)%NAME%OBJ=myName

       if_test: IF (TRIM(MyChannel) /= 'Test')   THEN   ! Test array
          CALL get_channel_object (status, TRIM(MyChannel), TRIM(myName)&
               , p4=CL(ic)%CPLDATA(ii)%ptr)
          IF (status /=0) write(iouerr,*) 'MMD2WAY_PARENT ERROR ' &
               , TRIM(MyChannel),' ',TRIM(myName)
          CALL channel_halt(substr, status)

          ! make sure that parent fields are present in restart files
          ! for restart frequency independent results
          CALL set_channel_object_restreq(status, TRIM(MyChannel), TRIM(myName))
          CALL channel_halt(substr, status)

          ! DEFINE AXIS STRING
          ! 1. get representation id
          CALL get_channel_object_info(status               &
               , TRIM(MyChannel),TRIM(myName), reprid=reprid)
          CALL channel_halt(substr, status)
          ! 2. get axis string
          CALL get_representation_info(status, ' ', reprid &
               , name=repr_name, gdimlen=repr_gdimlen      &
               , axis=CL(ic)%CPLDATA(ii)%AXIS              &
               , ldimlen=CL(ic)%CPLDATA(ii)%ldimlen)
          CALL channel_halt(substr, status)

          ! SEnd unit if requested
          IF (l_SendUnit) THEN
             CALL get_attribute(status,TRIM(MyChannel), TRIM(myName) &
                  , 'units', c=unit)
             IF (status /= 0) unit = ' '
             CALL MMD_Inter_Bcast(unit, .TRUE., Child_Id)
          END IF

          ! extract Representation Info if necessary
          IF (TRIM(ADJUSTL(clntrepr)) == '#UNKNOWN') THEN
             CALL get_channel_object_dimvalue(status  &
                  ,  TRIM(MyChannel), TRIM(myName)    &
                  , data=levs, axis ='Z', levtype=levtype)
             CALL channel_halt(substr, status)
             ch_att = ''
             IF (TRIM(levtype) == 'heights') THEN
                ch_att =str(levs(1))
                DO ix = 2,SIZE(levs)
                   ch_att = TRIM(ch_att)//','//str(levs(ix))
                END DO
             END IF

             CALL MMD_P_Send_Repr(CL(ic)%CPLDATA(ii)%AXIS &
                  , repr_gdimlen, repr_name, ch_att, Child_Id)
             IF (ASSOCIATED(levs)) THEN
                DEALLOCATE(levs)
                NULLIFY(levs)
             END IF
          ENDIF
       ELSE ! if_test
          CALL MMD_testP_GetTestPtr(child_id, CL(ic)%CPLDATA(ii)%ptr&
               , CL(ic)%CPLDATA(ii)%AXIS, CL(ic)%CPLDATA(ii)%ldimlen)
       ENDIF if_test

       IF (ASSOCIATED(CL(ic)%CPLDATA(ii)%ptr))   THEN
          CALL MMD_P_Set_DataArray (ic, status, CL(ic)%CPLDATA(ii)%ldimlen &
               , CL(ic)%CPLDATA(ii)%AXIS, CL(ic)%CPLDATA(ii)%ptr)
       ELSE
          CALL error_bi('P4 CPLDATA POINTER MUST BE ASSOCIATED!', substr)
       END IF
    END DO

    RETURN

  END SUBROUTINE Define_data_arrays
!-------------------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                PARENT COUPLING ROUTINES
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ----------------------------------------------------------------------
  SUBROUTINE mmd2way_parent_read_nml_cpl_child(status, iou, NML_TAG, lgrh &
       , i_rmy_px, rdefpc, lfreeslip, lcpl_gs, itype_fw, icosexp, damprel &
       , itype_VI, rcf, rcf_in, ldiagonly)

    ! MESSy/SMCL
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check &
                                      , read_nml_close
    USE messy_mmd2way,            ONLY: modstr

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(OUT)         :: status     ! error status
    INTEGER, INTENT(IN)          :: iou        ! I/O unit
    CHARACTER(LEN=3), INTENT(IN) :: NML_TAG
    LOGICAL, INTENT(OUT)         :: lgrh       ! humidity coupling
    INTEGER, INTENT(OUT)         :: i_rmy_px   ! exclude rmy+x from coupling
    REAL(dp), INTENT(OUT)        :: rdefpc     !
    LOGICAL,  INTENT(OUT)        :: lfreeslip  ! unforced boundary layer
    LOGICAL,  INTENT(OUT)        :: lcpl_gs    ! CH2P coupling in global start
    INTEGER,  INTENT(OUT)        :: itype_fw   ! weight function type
    INTEGER,  INTENT(OUT)        :: icosexp    ! exponent for itype_fw = 1
    REAL(dp), INTENT(OUT)        :: damprel    ! length of interface region
                                               !  itype_fw == 2
    INTEGER,  INTENT(OUT)        :: itype_VI   ! vertical integration type
    INTEGER,  INTENT(OUT)        :: RCF
    INTEGER,  INTENT(OUT)        :: RCF_IN
    LOGICAL,  INTENT(OUT)        :: ldiagonly

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'parent_read_nml_cpl_child'

    ! LOCAL
    LOGICAL          :: lex      ! file exists ?
    INTEGER          :: fstat    ! file status
    INTEGER          :: ii
    CHARACTER(LEN=3) :: INSTANCE = '   '

    NAMELIST /CPL_PAR_CHILD/ INSTANCE, lgrh, i_rmy_px, PFIELD &
         , rdefpc, lfreeslip, lcpl_gs, itype_fw, icosexp, damprel, itype_VI &
         , RCF, RCF_IN, ldiagonly

    status    = 1
    ! PRE-SET NAMELIST
    lgrh      = .FALSE.
    i_rmy_px  = 0
    rdefpc    = 30000._dp
    lfreeslip = .FALSE.
    lcpl_gs   = .FALSE.
    itype_fw  = 2
    icosexp   = 14
    damprel   = 0.2_dp
    itype_VI  = 1
    RCF       = 10000
    RCF_IN    = 10000
    ldiagonly = .FALSE.

    CALL read_nml_open(lex, substr, iou, 'CPL_PAR_CHILD', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist


    DO ! 64 = max number of coupling instances

       DO ii=1,NMAX_P_EXCH
          PFIELD(ii)%Parentname%CHA     = ' '
          PFIELD(ii)%Parentname%OBJ     = ' '
          PFIELD(ii)%Parenttend%CHA     = ' '
          PFIELD(ii)%Parenttend%OBJ     = ' '
          PFIELD(ii)%Childname%CHA      = ' '
          PFIELD(ii)%Childname%OBJ      = ' '
       END DO

       READ(iou, NML=CPL_PAR_CHILD, IOSTAT=fstat)

       CALL read_nml_check(fstat, substr, iou, 'CPL_PAR_CHILD', modstr)
       IF (fstat /= 0)  THEN
          IF (NML_TAG == 'DEF') THEN
             status = fstat
          ELSE
             status = 444 ! no namelist available or error while reading namelist
          ENDIF
          RETURN
       END IF
       ! in case of first call read only once, as first namelist
       ! is used as default namelist
       IF (NML_TAG == 'DEF') EXIT

       ! CORRECT NAMELIST FOUND
       IF (INSTANCE == NML_TAG) EXIT

    ENDDO

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES

    WRITE(*,*) NML_TAG, '....................................................'
    WRITE(*,*) NML_TAG, ' USE GENERALISED HUMIDITY? ', lgrh
    WRITE(*,*) NML_TAG, '....................................................'

    WRITE(*,*) NML_TAG, ' EXCLUDE RMY FROM COUPLING? i_rmy_px = ', i_rmy_px
    WRITE(*,*) NML_TAG, '....................................................'

    WRITE(*,*) NML_TAG, '....................................................'
    WRITE(*,*) NML_TAG, ' USE PCONTROL_FI ', rdefpc
    WRITE(*,*) NML_TAG, '....................................................'

    WRITE(*,*) NML_TAG, '....................................................'
    WRITE(*,*) NML_TAG, ' FREE SLIP in BOUNDARY LAYER? ', lfreeslip
    WRITE(*,*) NML_TAG, '....................................................'

    WRITE(*,*) NML_TAG, '....................................................'
    WRITE(*,*) NML_TAG, ' CHILD-PARENT COUPLING in GLOBAL START? ', lcpl_gs
    WRITE(*,*) NML_TAG, '....................................................'

    WRITE(*,*) NML_TAG, '....................................................'
    WRITE(*,*) NML_TAG, ' use weight function  ', itype_fw
    WRITE(*,*) NML_TAG, '....................................................'
    IF (itype_fw == 1) THEN
       WRITE(*,*) NML_TAG, '...................................................'
       WRITE(*,*) NML_TAG, ' use COSINUS exponent ',  icosexp
       WRITE(*,*) NML_TAG, '...................................................'
    ELSE IF (itype_fw == 2) THEN
       WRITE(*,*) NML_TAG, '...................................................'
       WRITE(*,*) NML_TAG, ' SIZE of damping zone ', damprel
       WRITE(*,*) NML_TAG, '...................................................'
    END IF

    WRITE(*,*) NML_TAG, '....................................................'
    WRITE(*,*) NML_TAG, ' use vertical integration (C2E) ', itype_VI
    WRITE(*,*) NML_TAG, '....................................................'

    WRITE(*,*) NML_TAG, '....................................................'
    WRITE(*,*) NML_TAG, ' INTEGER scaling factor child  grid ', RCF
    WRITE(*,*) NML_TAG, ' INTEGER scaling factor parent grid ', RCF_IN
    WRITE(*,*) NML_TAG, '....................................................'

    IF (ldiagonly) THEN
       WRITE(*,*) NML_TAG,'....................................................'
       WRITE(*,*) NML_TAG,'       ONLY DIAGNOSTIC C2E EXCHANGE REQUIRED '
       WRITE(*,*) NML_TAG,'....................................................'
    ENDIF

    cplfield_loop: DO ii = 1, NMAX_P_EXCH

       IF (TRIM(PFIELD(ii)%Parentname%CHA) == '') CYCLE
       IF (TRIM(PFIELD(ii)%Parentname%OBJ) == '') CYCLE
       IF (TRIM(PFIELD(ii)%Childname%CHA)  == '') CYCLE
       IF (TRIM(PFIELD(ii)%Childname%OBJ)  == '') CYCLE

       WRITE(*,*) NML_TAG, '....................................................'
       WRITE(*,*) NML_TAG, ' NUMBER OF CPL-FIELD: ',ii
       WRITE(*,*) NML_TAG, ' CHANNEL PARENT         = ' &
            , TRIM(PFIELD(ii)%Parentname%CHA)
       WRITE(*,*) NML_TAG, ' CHANNEL OBJECT PARENT  = ' &
            , TRIM(PFIELD(ii)%Parentname%OBJ)
       IF (TRIM(PFIELD(ii)%Parenttend%OBJ) /= '') THEN
          WRITE(*,*) NML_TAG, ' CHANNEL PARENT         = ' &
               , TRIM(PFIELD(ii)%Parenttend%CHA)
          WRITE(*,*) NML_TAG, ' CHANNEL OBJECT PARENT  = ' &
               , TRIM(PFIELD(ii)%Parenttend%OBJ)
       END IF

       WRITE(*,*) NML_TAG, ' CHANNEL CHILD         = ' &
            , TRIM(PFIELD(ii)%Childname%CHA)
       WRITE(*,*) NML_TAG, ' CHANNEL OBJECT CHILD  = ' &
            , TRIM(PFIELD(ii)%Childname%OBJ)

       WRITE(*,*) NML_TAG, ' REPRESENTATION         = ',  PFIELD(ii)%REPR
       WRITE(*,*) NML_TAG, ' INTERPOLATION METHOD   = ',  PFIELD(ii)%interp
       WRITE(*,*) NML_TAG, ' APPLICATION factor     = ',  PFIELD(ii)%fac
       WRITE(*,*) NML_TAG, ' APPLICATION METHOD     = ',  PFIELD(ii)%appl
       WRITE(*,*) NML_TAG, '....................................................'
    END DO cplfield_loop

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE mmd2way_parent_read_nml_cpl_child
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE interpret_parent_namelist

    ! MESSy/BMIL
    USE messy_main_mpi_bi,           ONLY: p_io, p_bcast, p_parallel_io
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! MESSy/SMCL
    USE messy_main_tools,      ONLY: match_wild
    USE messy_main_channel,    ONLY: get_channel_info

    IMPLICIT NONE

    INTRINSIC :: INDEX, SIZE, TRIM

    CHARACTER(LEN=*), PARAMETER            :: substr= 'interpret_parent_namelist'
    INTEGER, PARAMETER                     :: LMAX_EXCH = 10000
    TYPE(T_EXCH_P_IO), DIMENSION(LMAX_EXCH):: LFIELD
    LOGICAL                                :: l_wild_exist = .FALSE.
    CHARACTER(LEN=STRLEN_OBJECT), DIMENSION(:), POINTER  :: objname => NULL()
    INTEGER :: ij       ! loop over object list in case of wild card
    INTEGER :: ii       ! loop over exchange fields per child
    INTEGER :: ic       ! loop over childs
    INTEGER :: status   ! error status
    INTEGER :: NEXCH    ! counter number of valid exchange fields

    status = 1

    child_loop_02: DO ic = 1, NumberOfChildren

       IF (p_parallel_io) THEN

          ! COUNT  COUPLE FIELDS
          NEXCH = 0

          DO ii=1, NMAX_P_EXCH
             IF (TRIM(PIO(ic)%FIELD(ii)%Parentname%CHA) == '') CYCLE
             IF (TRIM(PIO(ic)%FIELD(ii)%Parentname%OBJ) == '') CYCLE
             IF (TRIM(PIO(ic)%FIELD(ii)%Childname%CHA)  == '') CYCLE
             IF (TRIM(PIO(ic)%FIELD(ii)%Childname%OBJ)  == '') CYCLE
             !

             l_wild_exist = (INDEX(PIO(ic)%FIELD(ii)%Parentname%OBJ, '?') > 0) &
                  .OR. (INDEX(PIO(ic)%FIELD(ii)%Parentname%OBJ, '*') > 0)
             IF (l_wild_exist) THEN
                ! WILDCARD NOT POSSIBLE FOR THE MMD2WAY_PARENT CHANNEL ITSELF
                IF (TRIM(PIO(ic)%FIELD(ii)%Parentname%CHA) == 'mmd2way_parent')&
                     CALL error_bi(substr, &
                     'wildcard for parent modules not allowed')

                CALL get_channel_info(status &
                     , TRIM(PIO(ic)%FIELD(ii)%Parentname%CHA), ONAMES=objname)
                CALL channel_halt(substr,status)
                DO ij = 1, SIZE(objname)
                   IF ( match_wild(TRIM(PIO(ic)%FIELD(ii)%Parentname%OBJ) &
                        ,TRIM(objname(ij)))) THEN
                      NEXCH   = NEXCH + 1
                      ! FIRST COPY SETTINGS
                      LFIELD(NEXCH) = PIO(ic)%FIELD(ii)
                      ! RESET WITH OBJECT SPECIFIC STUFF
                      LFIELD(NEXCH)%Parentname%OBJ = objname(ij)
                      LFIELD(NEXCH)%Childname%OBJ  = objname(ij)
                   ENDIF
                ENDDO
             ELSE
                NEXCH = NEXCH + 1
                LFIELD(NEXCH) = PIO(ic)%FIELD(ii)
             ENDIF
          ENDDO

          ! ALLOCATE CPLDATA STRUCT TO ACTUAL NUMBER OF COUPLING FIELDS
          ALLOCATE (PAR(ic)%CPLDATA(NEXCH))
          PAR(ic)%NEXCH = NEXCH

          WRITE(*,*) '.........................................................'
          WRITE(*,*) ' NUMBER OF Parent COUPLING-FIELDS: ',NEXCH
          WRITE(*,*) '.........................................................'

          DO ii=1, NEXCH
             PAR(ic)%CPLDATA(ii)%Parentname%CHA = LFIELD(ii)%Parentname%CHA
             PAR(ic)%CPLDATA(ii)%Parentname%OBJ = LFIELD(ii)%Parentname%OBJ
             PAR(ic)%CPLDATA(ii)%Parenttend%CHA = LFIELD(ii)%Parenttend%CHA
             PAR(ic)%CPLDATA(ii)%Parenttend%OBJ = LFIELD(ii)%Parenttend%OBJ
             PAR(ic)%CPLDATA(ii)%Childname%CHA  = LFIELD(ii)%Childname%CHA
             PAR(ic)%CPLDATA(ii)%Childname%OBJ  = LFIELD(ii)%Childname%OBJ
             PAR(ic)%CPLDATA(ii)%REPR   = LFIELD(ii)%REPR
             PAR(ic)%CPLDATA(ii)%INTERP = LFIELD(ii)%INTERP
             PAR(ic)%CPLDATA(ii)%APPL   = LFIELD(ii)%APPL
             PAR(ic)%CPLDATA(ii)%FAC    = LFIELD(ii)%FAC
             WRITE(*,*) '......................................................'
             WRITE(*,*) ' NUMBER OF EXCHANGE-FIELD: ',ii
             WRITE(*,*) ' CHANNEL - OBJECT PARENT = '                &
                  , TRIM(PAR(ic)%CPLDATA(ii)%Parentname%CHA),'  -  ' &
                  , TRIM(PAR(ic)%CPLDATA(ii)%Parentname%OBJ)
             IF (TRIM(PAR(ic)%CPLDATA(ii)%Parenttend%OBJ) /= '') THEN
                WRITE(*,*) ' CHANNEL - OBJECT PARENT TENDENCY = '       &
                     , TRIM(PAR(ic)%CPLDATA(ii)%Parenttend%CHA),'  -  ' &
                     , TRIM(PAR(ic)%CPLDATA(ii)%Parenttend%OBJ)
             END IF
             WRITE(*,*) ' CHANNEL - OBJECT CHILD = ',              &
                  TRIM(PAR(ic)%CPLDATA(ii)%Childname%CHA), '  -  ' &
                  ,TRIM(PAR(ic)%CPLDATA(ii)%Childname%OBJ)
             WRITE(*,*) ' REPRESENTATION       = ',  PAR(ic)%CPLDATA(ii)%REPR
             WRITE(*,*) ' INTERPOLATION METHOD = ',  PAR(ic)%CPLDATA(ii)%INTERP
             WRITE(*,*) ' APPLICATION METHOD   = ',  PAR(ic)%CPLDATA(ii)%APPL
             WRITE(*,*) ' APPLICATION FACTOR   = ',  PAR(ic)%CPLDATA(ii)%FAC
             WRITE(*,*) '......................................................'
          ENDDO

          PAR(ic)%lgrhum = PIO(ic)%lgrhum
          WRITE(*,*) '.........................................................'
          WRITE(*,*) ' USE GENERALISED HUMIDITY? ',PAR(ic)%lgrhum
          WRITE(*,*) '.........................................................'
          PAR(ic)%i_rmy_px = PIO(ic)%i_rmy_px
          WRITE(*,*) '.........................................................'
          WRITE(*,*) ' EXCLUDE RMY FROM COUPLING? i_rmy_px = ', PAR(ic)%i_rmy_px
          WRITE(*,*) '.........................................................'
          PAR(ic)%pcontrol_fi = PIO(ic)%pcontrol_fi
          WRITE(*,*) '.........................................................'
          WRITE(*,*) ' CONTROL GEOPOTENTIAL PCONTROL_FI ', PAR(ic)%pcontrol_fi
          WRITE(*,*) '.........................................................'

          PAR(ic)%lfreeslip = PIO(ic)%lfreeslip
          WRITE(*,*) '.........................................................'
          WRITE(*,*) ' FREE SLIP in BOUNDARY LAYER? ', PAR(ic)%lfreeslip
          WRITE(*,*) '.........................................................'

          PAR(ic)%lcpl_gs = PIO(ic)%lcpl_gs
          WRITE(*,*) '.........................................................'
          WRITE(*,*) ' CHILD-PARENT COUPLING in GLOBAL_START? ' &
               , PAR(ic)%lcpl_gs
          WRITE(*,*) '.........................................................'

          PAR(ic)%itype_fw = PIO(ic)%itype_fw
          WRITE(*,*) '.........................................................'
          WRITE(*,*) ' USE weight function ', PAR(ic)%itype_fw
          WRITE(*,*) '.........................................................'
          PAR(ic)%icosexp = PIO(ic)%icosexp
          PAR(ic)%damprel = PIO(ic)%damprel
          IF ( PAR(ic)%itype_fw == 1) THEN
             WRITE(*,*) '......................................................'
             WRITE(*,*) ' Use COSINUS exponent of ', PAR(ic)%icosexp
             WRITE(*,*) '......................................................'
          ELSE IF ( PAR(ic)%itype_fw == 2) THEN
             WRITE(*,*) '......................................................'
             WRITE(*,*) ' relative width of damping zone ', PAR(ic)%damprel
             WRITE(*,*) '......................................................'
          ELSE IF ( PAR(ic)%itype_fw == 0) THEN
             WRITE(*,*) '......................................................'
             WRITE(*,*) ' set weight function 1 everywhere '
             WRITE(*,*) '......................................................'
          END IF

          PAR(ic)%itype_VI = PIO(ic)%itype_VI
          WRITE(*,*) '.........................................................'
          WRITE(*,*) ' USE vertical integration ', PAR(ic)%itype_VI
          WRITE(*,*) '.........................................................'

          PAR(ic)%RCF = PIO(ic)%RCF
          WRITE(*,*) '.........................................................'
          WRITE(*,*) ' USE scaling factor CHILD GRID ', PAR(ic)%RCF
          WRITE(*,*) '.........................................................'

          PAR(ic)%RCF_IN = PIO(ic)%RCF_IN
          WRITE(*,*) '.........................................................'
          WRITE(*,*) ' USE scaling factor PARENT GRID ', PAR(ic)%RCF_IN
          WRITE(*,*) '.........................................................'

          PAR(ic)%ldiagonly = PIO(ic)%ldiagonly
          WRITE(*,*) '.........................................................'
          WRITE(*,*) ' SOLELY EXCHANGE OF DIAGNOSTIC FIELDS',PAR(ic)%ldiagonly
          WRITE(*,*) '.........................................................'

       ENDIF

       ! DISTRIBUTE VALUES
       CALL p_bcast(NEXCH, p_io)
       IF (.NOT. p_parallel_io) THEN
          ALLOCATE (PAR(ic)%CPLDATA(NEXCH))
          PAR(ic)%NEXCH = NEXCH
       ENDIF

       DO  ii=1,NEXCH
          CALL p_bcast(PAR(ic)%CPLDATA(ii)%Parentname%CHA, p_io)
          CALL p_bcast(PAR(ic)%CPLDATA(ii)%Parentname%OBJ, p_io)
          CALL p_bcast(PAR(ic)%CPLDATA(ii)%Parenttend%CHA, p_io)
          CALL p_bcast(PAR(ic)%CPLDATA(ii)%Parenttend%OBJ, p_io)
          CALL p_bcast(PAR(ic)%CPLDATA(ii)%Childname%CHA,  p_io)
          CALL p_bcast(PAR(ic)%CPLDATA(ii)%Childname%OBJ,  p_io)
          CALL p_bcast(PAR(ic)%CPLDATA(ii)%REPR,           p_io)
          CALL p_bcast(PAR(ic)%CPLDATA(ii)%INTERP,         p_io)
          CALL p_bcast(PAR(ic)%CPLDATA(ii)%APPL,           p_io)
          CALL p_bcast(PAR(ic)%CPLDATA(ii)%FAC,            p_io)
          IF (PAR(ic)%CPLDATA(ii)%INTERP /= 1)                                 &
             CALL error_bi('only conservative interpolation implemented so far'&
                  , substr)
          ! treatment of input fields (appl=0)
          IF (PAR(ic)%CPLDATA(ii)%appl<0 .OR. PAR(ic)%CPLDATA(ii)%appl>2) &
             CALL error_bi('invalid application method', substr)
          IF (PAR(ic)%CPLDATA(ii)%appl==0) THEN
             IF (PAR(ic)%CPLDATA(ii)%Parentname%CHA(1:7) == 'mmd2way') THEN
                IF (p_parallel_io) write (*,*) ' mask field: '       &
                     , TRIM(PAR(ic)%CPLDATA(ii)%Parentname%CHA),'  ' &
                     , TRIM(PAR(ic)%CPLDATA(ii)%Parentname%OBJ)
             ELSE
                IF (p_parallel_io) write (*,*) ' masking field '     &
                     , TRIM(PAR(ic)%CPLDATA(ii)%Parentname%CHA),'  ' &
                     , TRIM(PAR(ic)%CPLDATA(ii)%Parentname%OBJ)      &
                     , ' NOT POSSIBLE!'
                CALL error_bi('masking of object not possible', substr)
             ENDIF
          ENDIF
       ENDDO

    ENDDO child_loop_02

    DEALLOCATE(PIO)

    status = 0

  END SUBROUTINE interpret_parent_namelist
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE Setup_ParData_Exchange(child_id,ic)

    ! MESSy/BMIL
    USE messy_main_transform_bi,  ONLY: locate_in_decomp
    USE messy_main_mpi_bi,        ONLY: p_pe
    ! MMD
    USE mmd_parent,               ONLY: MMD_Recv_from_Child     &
                                      , MMD_Send_to_Child       &
                                      , MMD_P_Get_ParIndexList

    IMPLICIT NONE

    INTRINSIC :: SIZE

    INTEGER,INTENT(IN)    :: ic
    INTEGER,INTENT(IN)    :: child_id

    ! LOCAL
    INTEGER               :: i,j,k, status
    INTEGER               :: icount           ! Index cout
    INTEGER               :: icnt(1)          ! summy for MPI communication
    INTEGER, DIMENSION(3) :: size_of_array
    INTEGER, DIMENSION(3) :: cgrid
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: index_list

#ifdef COSMO
    CHARACTER(LEN=*), PARAMETER :: substr='Setup_ParData_Exchange_cosmo'
#endif
#ifdef ECHAM5
    CHARACTER(LEN=*), PARAMETER :: substr='Setup_ParData_Exchange_echam'
#endif

    IF (PAR(ic)%NEXCH == 0) RETURN
    if_pe: IF (p_pe == 0) THEN
      CALL MMD_Recv_from_Child (child_id, size_of_array, 3, 0, 140, status)
      IF (status /= 0) CALL error_bi('MMD_Recv_from_Child ERROR 01', substr)

      ALLOCATE(PAR(ic)%all_plon(size_of_array(1) &
           ,size_of_array(2),size_of_array(3)))
      ALLOCATE(PAR(ic)%all_plat(size_of_array(1) &
           ,size_of_array(2),size_of_array(3)))

      CALL MMD_Recv_from_Child (child_id, PAR(ic)%all_plon &
           , SIZE(PAR(ic)%all_plon), 0, 141, status)
      IF (status /= 0) CALL error_bi('MMD_Recv_from_Child ERROR 02', substr)
      CALL MMD_Recv_from_Child (child_id, PAR(ic)%all_plat &
           , SIZE(PAR(ic)%all_plat), 0, 142, status)
      IF (status /= 0) CALL error_bi('MMD_Recv_from_Child ERROR 03', substr)


      ALLOCATE (index_list(6, SIZE(PAR(ic)%all_plon)))

      ! Create list of required PARENT grid boxes
      icount = 0
      DO k=1,size_of_array(3)
         DO j=1,size_of_array(2)
            DO i=1,size_of_array(1)
               ! Legal Value
               if_legal: IF (PAR(ic)%all_plon(i,j,k) > -1000._dp .AND. &
                    PAR(ic)%all_plat(i,j,k) > -1000._dp) THEN

                  CALL locate_in_decomp (status, PAR(ic)%all_plon(i,j,k)     &
                       , PAR(ic)%all_plat(i,j,k),cgrid(1),cgrid(2),cgrid(3)  &
#ifdef COSMO
                       , linclbnd=.TRUE., lprinterr=.TRUE., lrot=.TRUE.      &
#endif
                      )


                  IF (status /= 0) THEN
                     write(iouerr,*) 'POINT NOT LOCATED IN PARENT DOMAIN '  &
                          , i, j,icount+1,' cgrid ',cgrid, ' CHILD_PE ',k-1 &
                          , PAR(ic)%all_plon(i,j,k), PAR(ic)%all_plat(i,j,k )
                     CALL error_bi(substr, 'POINT NOT LOCATED IN PARENT DOMAIN')
                  ENDIF
                  icount = icount+1
                  index_list(1,icount) = cgrid(2) ! First index in Parent Array
                  index_list(2,icount) = cgrid(3) ! Second index in Parent Array
                  index_list(3,icount) = i        ! Lon Index child coarse grid
                  index_list(4,icount) = j        ! Lat Index child coarse grid
                  index_list(5,icount) = k-1      ! PE Number child
                  index_list(6,icount) = cgrid(1) ! PE Number parent

                  IF (index_list(6,icount) == -1 )                    &
                       write(iouerr,*) 'ERROR ParIndList ', i,j,k     &
                       , index_list(:,icount),PAR(ic)%all_plon(i,j,k) &
                       , PAR(ic)%all_plat(i,j,k)
               END IF if_legal
            END DO
         END DO
      END DO

      icnt(1) = icount
      CALL MMD_Send_to_Child (child_id, icnt(1:1), 1, 0, 144, status)
      IF (status /= 0) CALL error_bi('MMD_Send_to_Child 01', substr)
      CALL MMD_Send_to_Child (child_id, index_list(:,1:icount) &
                                        , 6*icount, 0, 145, status)
      IF (status /= 0) CALL error_bi('MMD_Send_to_Child 02', substr)

      DEALLOCATE(index_list)
   END IF if_pe

   CALL MMD_P_Get_ParIndexList(ic, PAR(ic)%frac, Par(ic)%wf)

   DO i=1,SIZE(PAR(ic)%frac,1)
      DO j=1,SIZE(PAR(ic)%frac,2)
         IF (PAR(ic)%frac(i,j) > 1.005_dp) THEN
            write(iouerr,*) 'PARENT PE ', p_pe &
                 ,'BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB'
            write(iouerr,*) 'PARENT PE ', p_pe , ' weights larger' &
                                    , i,j, PAR(ic)%frac(i,j)
            write(iouerr,*) 'PARENT PE ', p_pe    &
                 ,'BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB'
         ENDIF
      END DO
   END DO

   ! DEFINE MASK (allows for multiplication instead of if-statements
   PAR(ic)%mask(:,:) = 0._dp
   DO i=1,SIZE(PAR(ic)%frac,1)
      DO j=1,SIZE(PAR(ic)%frac,2)
         IF (PAR(ic)%frac(i,j) > 0.99999_dp) THEN
            PAR(ic)%mask(i,j) = 1._dp
         ENDIF
      END DO
   END DO

   RETURN
 END SUBROUTINE Setup_ParData_Exchange
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE  mmd2way_parent_couple_gp(ic, flag)

    ! MESSy/BMIL
    USE messy_main_grid_def_mem_bi, ONLY: ngpblks, nproma
#ifdef ECHAM5
    USE messy_main_grid_def_mem_bi, ONLY: npromz
#endif
#ifdef COSMO
    USE messy_main_grid_def_mem_bi,  ONLY: jstartpar, jendpar &
                                         , nlev, ie, je
    USE messy_main_grid_def_bi,      ONLY: lperi_x, lperi_y, l2dim
    USE messy_main_data_bi,          ONLY: nnew

    USE messy_main_mpi_bi,           ONLY: icomm_cart, sendbuf, isendbuflen  &
                                         , imp_reals, nboundlines            &
                                         , my_cart_neigh, ncomm_type         &
                                         , ldatatypes, num_compute           &
                                         , exchg_boundaries
#endif
    USE messy_main_tracer_mem_bi, ONLY: xtte, xtm1
    ! MESSy/SMCL
    USE messy_main_timer,         ONLY: time_step_len

    IMPLICIT NONE

    ! INPUT IO
    INTEGER, INTENT(IN) :: ic  ! child ID
    INTEGER, INTENT(IN) :: flag
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr ='mmd2way_parent_couple_gp'
    INTEGER  :: ii, jrow, jk  ! loop indices
    INTEGER  :: jkkm          ! loop index for reduced height coupling
    INTEGER  :: kmin, kmax    ! indices vertical extent
    INTEGER  :: idt           ! local tracer index
    INTEGER  :: tlev          ! time level, required for COSMO
    INTEGER  :: kproma, jp
    REAL(dp) ::  tend, fc
#ifdef COSMO
    ! required for exchg_boundaries
    INTEGER  :: kdims(24), status, size3, icount
    CHARACTER(LEN=80) :: yerrmsg = ' '

    icount = 1111
#endif

    array_loop: DO ii = 1, Par(ic)%NEXCH
       IF (TRIM(Par(ic)%CPLData(ii)%parentname%obj) == 'Test_Ar') CYCLE
       IF (Par(ic)%CPLData(ii)%appl > 1 )                         CYCLE

       IF (.NOT. Par(ic)%CPLData(ii)%fac > 0._dp)  CYCLE
       idt = Par(ic)%CPLData(ii)%idt

       ! CHANGE TRACER / TENDENCIES before final integration
       IF (Par(ic)%CPLData(ii)%lte .AND. flag /= 1) CYCLE
       ! CHANGE FULL FIELDS after FINAL INTEGRATION
       IF (.NOT. Par(ic)%CPLData(ii)%lte .AND. flag == 1) CYCLE

#ifdef COSMO
       IF (Par(ic)%CPLData(ii)%ltimedep) THEN
          tlev = nnew
       ELSE
          tlev = 1
       END IF
#else
       tlev = 1
#endif
       if_trac: IF (idt > 0) THEN
          ! array  is tracer
          ! enable coupling in reduced height range
          kmin = Par(ic)%kmin
          kmax = SIZE(Par(ic)%CPLData(ii)%target,_IZ_XYZ__)

          jrow1: DO jrow = 1, ngpblks
             jk1: DO jk = kmin + 1,kmax
                jkkm = jk - kmin + 1
                kproma = nproma
#ifdef ECHAM5
                IF (jrow == ngpblks)  kproma = npromz
#endif
                jp1: DO jp = 1, kproma
                   if_frac: IF (Par(ic)%frac(jp,jrow) > 0.99999_dp) THEN
                  ! absolute value of tracer mixing ratio was exchanged
                  !
                  ! calculate tendency:
                  ! X = (1-fc) Xe + f Xc
                  ! Xm1 + dX dt = (1-fc) (Xm1 + dXe dt)  + fc Xc
                  !       dX dt = (1-fc) dXe dt - fc Xm1 + fc Xc
                  !       dX    = (1-fc) dXe + fc (Xc -Xm1) / dt

                      fc = Par(ic)%wf(jp,jrow)               &
                           ! magnitude of nudging
                           * Par(ic)%CPLData(ii)%fac         &
                           ! coupling in reduced height range
                           * Par(ic)%kminfac(jk)

                      tend = (Par(ic)%CPLData(ii)%ptr(_RI_XYZ__(jp,jrow,jkkm),1) &
                              - xtm1(_RI_XYZN_(jp,jrow,jk,idt))) / time_step_len
                      Par(ic)%CPLData(ii)%target(_RI_XYZ__(jp,jrow,jk),1) = fc * tend
                      xtte(_RI_XYZN_(jp,jrow,jk,idt)) = (1._dp- fc) * xtte(_RI_XYZN_(jp,jrow,jk,idt)) &
                           + fc * tend
                  END IF if_frac
                END DO jp1
             END DO jk1
          END DO jrow1
#ifdef COSMO
          kdims(1:24) = &
               (/nlev,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
          CALL exchg_boundaries                                       &
               (100+icount,sendbuf,isendbuflen, imp_reals             &
               , icomm_cart, num_compute                              &
               , ie, je, kdims, jstartpar, jendpar                    &
               , nboundlines, nboundlines, my_cart_neigh              &
               , lperi_x, lperi_y, l2dim, 50000+icount                &
               , ldatatypes, ncomm_type, status, yerrmsg              &
               , xtte(_RI_XYZN_(:,:,:,idt)))
          IF (status /= 0) CALL error_bi(yerrmsg, substr)
          icount=icount+1
#endif
       ELSE ! if_trac
          ! ARRAY (not tracer)
          if_lte: IF (Par(ic)%CPLData(ii)%lte) THEN
             ! force tendency
             IF (Par(ic)%CPLData(ii)%rank == 2) THEN
                if_ps: IF (Par(ic)%CPLData(ii)%Parentname%OBJ(1:4)=='alps') THEN
                   kproma = nproma

                   jrow2a: DO jrow = 1, ngpblks
#ifdef ECHAM5
                      IF (jrow == ngpblks)  kproma = npromz
#endif
                      jp2a: DO jp = 1, kproma
                         if_frac1a: IF (Par(ic)%frac(jp,jrow) > 0.99999_dp) THEN

                            fc = Par(ic)%wf(jp,jrow)*Par(ic)%CPLData(ii)%fac
                            tend = &
                              (log(Par(ic)%CPLData(ii)%ptr(jp,jrow,1,1)) &
                              - Par(ic)%CPLData(ii)%targetm1(jp,jrow,1,1))  &
                              / time_step_len
                            Par(ic)%CPLData(ii)%target(jp,jrow,1,1) = &
                                 (1._dp- fc) &
                                    * Par(ic)%CPLData(ii)%target(jp,jrow,1,1) &
                                + fc * tend
                         ENDIF if_frac1a
                      END DO jp2a
                   END DO jrow2a
#ifdef COSMO
                   kdims(1:24) = &
                        (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
                   CALL exchg_boundaries                                       &
                        (100+icount,sendbuf,isendbuflen, imp_reals             &
                        , icomm_cart, num_compute                              &
                        , ie, je, kdims, jstartpar, jendpar                    &
                        , nboundlines, nboundlines, my_cart_neigh              &
                        , lperi_x, lperi_y, l2dim, 50000+icount                &
                        , ldatatypes, ncomm_type, status, yerrmsg              &
                        , Par(ic)%CPLData(ii)%target(:,:,:,1))
                   IF (status /= 0) CALL error_bi(yerrmsg, substr)
                   icount=icount+1
#endif
                ELSE ! if_ps => no ps
                   kproma = nproma

                   jrow2: DO jrow = 1, ngpblks
#ifdef ECHAM5
                      IF (jrow == ngpblks)  kproma = npromz
#endif
                      jp2: DO jp = 1, kproma
                         if_frac1: IF (Par(ic)%frac(jp,jrow) > 0.99999_dp) THEN
                            fc = Par(ic)%wf(jp,jrow)*Par(ic)%CPLData(ii)%fac
                            tend =                                             &
                              (Par(ic)%CPLData(ii)%ptr(jp,jrow,1,1)            &
                              - Par(ic)%CPLData(ii)%targetm1(jp,jrow,tlev,1))  &
                              / time_step_len
                            Par(ic)%CPLData(ii)%target(jp,jrow,1,1) =          &
                                 (1._dp- fc)                                   &
                                  * Par(ic)%CPLData(ii)%target(jp,jrow,tlev,1) &
                                + fc * tend
                         ENDIF if_frac1
                      END DO jp2
                   END DO jrow2
#ifdef COSMO
                   kdims(1:24) = &
                        (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
                   CALL exchg_boundaries                                       &
                        (100+icount,sendbuf,isendbuflen, imp_reals             &
                        , icomm_cart, num_compute                              &
                        , ie, je, kdims, jstartpar, jendpar                    &
                        , nboundlines, nboundlines, my_cart_neigh              &
                        , lperi_x, lperi_y, l2dim, 50000+icount                &
                        , ldatatypes, ncomm_type, status, yerrmsg              &
                        , Par(ic)%CPLData(ii)%target(:,:,tlev,1))
                   IF (status /= 0) CALL error_bi(yerrmsg, substr)
                   icount=icount+1
#endif
                END IF if_ps
             ELSE IF (Par(ic)%CPLData(ii)%rank == 3) THEN
                kmin = Par(ic)%kmin
                kmax = SIZE(Par(ic)%CPLData(ii)%target,_IZ_XYZ__)
                jrow3: DO jrow = 1, ngpblks
                   jk3: DO jk = kmin + 1,kmax
                   jkkm = jk - kmin + 1
                      kproma = nproma
#ifdef ECHAM5
                      IF (jrow == ngpblks)  kproma = npromz
#endif
                      jp3: DO jp = 1, kproma
                         if_frac1b: IF (Par(ic)%frac(jp,jrow) > 0.99999_dp) THEN
                           fc = Par(ic)%wf(jp,jrow)                &
                                 ! magnitude of nudging
                                 * Par(ic)%CPLData(ii)%fac         &
                                 ! coupling in reduced height range
                                 * Par(ic)%kminfac(jk)

                           tend = (Par(ic)%CPLData(ii)%ptr(_RI_XYZ__(jp,jrow,jkkm),1)     &
                               - Par(ic)%CPLData(ii)%targetm1(_RI_XYZ__(jp,jrow,jk),tlev))&
                               / time_step_len
                           Par(ic)%CPLData(ii)%target(_RI_XYZ__(jp,jrow,jk),tlev) =    &
                                 (1._dp- fc) *                              &
                                 Par(ic)%CPLData(ii)%target(_RI_XYZ__(jp,jrow,jk),tlev)&
                                       + fc  * tend
                         END IF if_frac1b
                      END DO jp3
                   END DO jk3
                END DO jrow3
#ifdef COSMO
                size3 = SIZE(Par(ic)%CPLData(ii)%target,_IZ_XYZ__)
                kdims(1:24) = &
                     (/size3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
                CALL exchg_boundaries                                       &
                     (100+icount,sendbuf,isendbuflen, imp_reals             &
                     , icomm_cart, num_compute                              &
                     , ie, je, kdims, jstartpar, jendpar                    &
                     , nboundlines, nboundlines, my_cart_neigh              &
                     , lperi_x, lperi_y, l2dim, 50000+icount                &
                     , ldatatypes, ncomm_type, status, yerrmsg              &
                     , Par(ic)%CPLData(ii)%target(:,:,:,1))
                IF (status /= 0) CALL error_bi(yerrmsg, substr)
                icount=icount+1
#endif
             ELSE
                CALL error_bi(substr, 'rank 4D not yet implemented')
             END IF

          ELSE ! if_lte : change full fields


            if_rank2b: IF (Par(ic)%CPLData(ii)%rank == 2) THEN
               ifappl1a: IF (Par(ic)%CPLData(ii)%appl == 1) THEN
                  kproma = nproma
                  jrow2b: DO jrow = 1, ngpblks
#ifdef ECHAM5
                     IF (jrow == ngpblks)  kproma = npromz
#endif
                     jp2b: DO jp = 1, kproma
                        if_frac2b: IF (Par(ic)%frac(jp,jrow) > 0.99999_dp) THEN

                           Par(ic)%CPLData(ii)%target(jp,jrow,tlev,1) =        &
                                Par(ic)%CPLData(ii)%target(jp,jrow,tlev,1)     &
                                * (1._dp -                                     &
                                Par(ic)%CPLData(ii)%fac * Par(ic)%wf(jp,jrow)) &
                                + Par(ic)%CPLData(ii)%ptr(jp,jrow,1,1)         &
                                ! magnitude of nudging
                                * Par(ic)%CPLData(ii)%fac * Par(ic)%wf(jp,jrow)
                        ENDIF if_frac2b
                     END DO jp2b
                  END DO jrow2b
#ifdef COSMO
                  kdims(1:24) = &
                       (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
                  CALL exchg_boundaries                                    &
                       (100+icount,sendbuf,isendbuflen, imp_reals          &
                       , icomm_cart, num_compute                           &
                       , ie, je, kdims, jstartpar, jendpar                 &
                       , nboundlines, nboundlines, my_cart_neigh           &
                       , lperi_x, lperi_y, l2dim, 50000+icount             &
                       , ldatatypes, ncomm_type, status, yerrmsg           &
                       , Par(ic)%CPLData(ii)%target(:,:,tlev,1))
                  IF (status /= 0) CALL error_bi(yerrmsg, substr)
                  icount=icount+1
#endif
               ! input fields appl=0 added
               ELSE
                  kproma = nproma
                  jrow5b: DO jrow = 1, ngpblks
#ifdef ECHAM5
                     IF (jrow == ngpblks)  kproma = npromz
#endif
                     jp5b: DO jp = 1, kproma
                        if_frac5b: IF (Par(ic)%frac(jp,jrow) > 0.99999_dp) THEN
                           Par(ic)%CPLData(ii)%target(jp,jrow,tlev,1) = &
                                Par(ic)%CPLData(ii)%ptr(jp,jrow,1,1)
                        ENDIF if_frac5b
                     END DO jp5b
                  END DO jrow5b
               ENDIF ifappl1a
#ifdef COSMO
                kdims(1:24) = &
                     (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
                CALL exchg_boundaries                                    &
                     (100+icount,sendbuf,isendbuflen, imp_reals          &
                     , icomm_cart, num_compute                           &
                     , ie, je, kdims, jstartpar, jendpar                 &
                     , nboundlines, nboundlines, my_cart_neigh           &
                     , lperi_x, lperi_y, l2dim, 50000+icount             &
                     , ldatatypes, ncomm_type, status, yerrmsg           &
                     , Par(ic)%CPLData(ii)%target(:,:,tlev,1))
                IF (status /= 0) CALL error_bi(yerrmsg, substr)
                icount=icount+1
#endif
            ELSE IF (Par(ic)%CPLData(ii)%rank == 3) THEN ! if_rank2
               ifappl1b: IF (Par(ic)%CPLData(ii)%appl == 1) THEN
                  kmin = Par(ic)%kmin
                  kmax = SIZE(Par(ic)%CPLData(ii)%target,_IZ_XYZ__)
                  jrow3b: DO jrow = 1, ngpblks
                     jk3b: DO jk = kmin + 1,kmax
                        jkkm = jk - kmin + 1
                        kproma = nproma
#ifdef ECHAM5
                        IF (jrow == ngpblks)  kproma = npromz
#endif
                        jp3b: DO jp = 1, kproma
                           IF (Par(ic)%frac(jp,jrow) > 0.99999_dp) THEN
                              Par(ic)%CPLData(ii)%target(_RI_XYZ__(jp,jrow,jk),tlev) =    &
                                   Par(ic)%CPLData(ii)%target(_RI_XYZ__(jp,jrow,jk),tlev) &
                                   * (1._dp - Par(ic)%CPLData(ii)%fac          &
                                                     * Par(ic)%kminfac(jk)     &
                                   * Par(ic)%wf(jp,jrow))                      &
                                   + Par(ic)%CPLData(ii)%ptr(_RI_XYZ__(jp,jrow,jkkm),1)   &
                                   ! magnitude of nudging
                                   * Par(ic)%CPLData(ii)%fac                   &
                                   ! coupling in reduced height range
                                   * Par(ic)%kminfac(jk)                       &
                                   ! weight function
                                   * Par(ic)%wf(jp,jrow)
                           END IF
                        END DO jp3b
                     END DO jk3b
                  END DO jrow3b
#ifdef COSMO
                size3 = SIZE(Par(ic)%CPLData(ii)%target,_IZ_XYZ__)
                kdims(1:24) = &
                     (/size3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
                CALL exchg_boundaries                                       &
                     (100+icount,sendbuf,isendbuflen, imp_reals             &
                     , icomm_cart, num_compute                              &
                     , ie, je, kdims, jstartpar, jendpar                    &
                     , nboundlines, nboundlines, my_cart_neigh              &
                     , lperi_x, lperi_y, l2dim, 50000+icount                &
                     , ldatatypes, ncomm_type, status, yerrmsg              &
                     , Par(ic)%CPLData(ii)%target(:,:,:,tlev))
                IF (status /= 0) CALL error_bi(yerrmsg, substr)
                icount=icount+1
#endif
               ELSE ! ifappl1b
                  kmin = Par(ic)%kmin
                  kmax = SIZE(Par(ic)%CPLData(ii)%target,_IZ_XYZ__)
                  jrow4b: DO jrow = 1, ngpblks
                     jk4b: DO jk = kmin + 1,kmax
                        jkkm = jk - kmin + 1
                        kproma = nproma
#ifdef ECHAM5
                        IF (jrow == ngpblks)  kproma = npromz
#endif
                        jp4b: DO jp = 1, kproma
                           IF (Par(ic)%frac(jp,jrow) > 0.99999_dp) THEN
                              Par(ic)%CPLData(ii)%target(_RI_XYZ__(jp,jrow,jk),tlev) =&
                                   Par(ic)%CPLData(ii)%ptr(_RI_XYZ__(jp,jrow,jkkm),1)
                           END IF
                        END DO jp4b
                     END DO jk4b
                  END DO jrow4b
               END IF ifappl1b
#ifdef COSMO
                size3 = SIZE(Par(ic)%CPLData(ii)%target,_IZ_XYZ__)
                kdims(1:24) = &
                     (/size3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
                CALL exchg_boundaries                                       &
                     (100+icount,sendbuf,isendbuflen, imp_reals     &
                     , icomm_cart, num_compute                              &
                     , ie, je, kdims, jstartpar, jendpar                    &
                     , nboundlines, nboundlines, my_cart_neigh              &
                     , lperi_x, lperi_y, l2dim, 50000+icount                &
                     , ldatatypes, ncomm_type, status, yerrmsg              &
                     , Par(ic)%CPLData(ii)%target(:,:,:,tlev))
                IF (status /= 0) CALL error_bi(yerrmsg, substr)
                icount=icount+1
#endif
             ELSE ! if_rank2
                CALL error_bi(substr, 'rank 4D not yet implemented')
             END IF if_rank2b
          ENDIF if_lte
       END IF if_trac ! tracer

    END DO array_loop

  END SUBROUTINE mmd2way_parent_couple_gp

  ! ----------------------------------------------------------------------
#else

  PUBLIC :: mmd2way_parent_initialize
  PUBLIC :: mmd2way_parent_init_memory
  PUBLIC :: mmd2way_parent_init_coupling
  PUBLIC :: mmd2way_parent_global_start
  PUBLIC :: mmd2way_parent_global_end
  PUBLIC :: mmd2way_parent_free_memory

CONTAINS

  SUBROUTINE mmd2way_parent_initialize
  END SUBROUTINE mmd2way_parent_initialize

  SUBROUTINE mmd2way_parent_init_memory
  END SUBROUTINE mmd2way_parent_init_memory

  SUBROUTINE mmd2way_parent_init_coupling
  END SUBROUTINE mmd2way_parent_init_coupling

  SUBROUTINE mmd2way_parent_global_start
  END SUBROUTINE mmd2way_parent_global_start

  SUBROUTINE mmd2way_parent_global_end(flag)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: flag
    IF (flag == 0) RETURN
  END SUBROUTINE mmd2way_parent_global_end

  SUBROUTINE mmd2way_parent_free_memory
  END SUBROUTINE mmd2way_parent_free_memory

#endif
END MODULE messy_mmd2way_parent_si
