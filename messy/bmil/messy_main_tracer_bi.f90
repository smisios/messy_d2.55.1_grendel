#include "messy_main_ppd_bi.inc"

! **********************************************************************+
MODULE messy_main_tracer_bi
! **********************************************************************+

  USE messy_main_tracer_mem_bi
  USE messy_main_blather_bi,      ONLY: start_message_bi, end_message_bi
  USE messy_main_constants_mem,   ONLY: FLAGGED_BAD
  USE messy_main_tracer_tools_bi, ONLY: tracer_halt
#ifdef ICON
  USE messy_main_tools,           ONLY: PTR_3D_ARRAY
  USE messy_main_tracer_tools_bi, ONLY: main_tracer_set_domain
#endif
  USE messy_main_tracer

  IMPLICIT NONE
  PRIVATE

  TYPE T_INIPVAR_IO
     CHARACTER(LEN=20)           :: basemodel = '' ! basemodel name
     CHARACTER(LEN=STRLEN_FNAME) :: pvar      = '' ! name of prog. var.
     CHARACTER(LEN=STRLEN_FNAME) :: inp       = '' ! name of input field
     INTEGER                     :: switch    = 0  ! switch for conversion
  END TYPE T_INIPVAR_IO
  !
  INTEGER, PARAMETER :: N_INI_MAX = 20
  !
  TYPE(T_INIPVAR_IO), DIMENSION(N_INI_MAX), SAVE :: ini_pvar

  ! GLOBAL NAMELIST SWITCHES (CPL)
  LOGICAL,  SAVE :: l_tracer_init = .TRUE.
  LOGICAL,  SAVE :: l_tracer_initfromrestart = .FALSE.
  ! only for LG->GP conversion (ATTILA)
  LOGICAL,  SAVE :: l_conv_lg2gp      = .TRUE.        ! do it?
  INTEGER,  SAVE :: i_conv_lg2gp_mode = 2             ! LG2GP_AVE
  REAL(DP), SAVE :: r_conv_lg2gp_fill = FLAGGED_BAD   ! fill with?
  LOGICAL,  SAVE :: l_conv_lg2gp_mc   = .FALSE.       ! mass conservarions?
  ! only for LG->GP conversion (CLaMS)
  LOGICAL,  SAVE :: l_conv_cl2gp      = .TRUE.        ! do it?
  INTEGER,  SAVE :: i_conv_cl2gp_mode = 2             ! CL2GP_AVE
  REAL(DP), SAVE :: r_conv_cl2gp_fill = FLAGGED_BAD   ! fill with?

#ifdef ICON
  TYPE(PTR_3D_ARRAY), DIMENSION(:),     POINTER :: grmass_conv    => NULL()
  TYPE(PTR_3D_ARRAY), DIMENSION(:),     POINTER :: grmassdry_conv => NULL()
#endif
  ! PUBLIC SUBROUTINES CALLED FROM BM(I)L
  ! ECHAM5/MESSy SPECIFIC
  PUBLIC :: main_tracer_setup         ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_initialize    ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_new_tracer    ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_init_memory   ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_init_coupling ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_init_tracer   ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_global_start  ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_beforeadv     ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_afteradv      ! CALLED FROM messy_main_control_e5/c4
  ! SPECIAL: CALLED BEFORE OUTPUT
  PUBLIC :: main_tracer_write_output  ! CALLED FROM messy_main_control_e5/c4
!!#D attila +
  !PRIVATE :: tracer_out_conv_lg2gp    ! CONVERT LG TRACERS TO GP FOR OUTPUT
!!#D attila -
!!#D clams +
  !PRIVATE :: tracer_out_conv_cl2gp    ! CONVERT CL TRACERS TO GP FOR OUTPUT
!!#D clams -
  !
  PUBLIC :: main_tracer_local_start   ! CALLED FROM messy_main_control_e5/c4
#if defined(ICON)
  PUBLIC :: main_tracer_turbdiff
  PUBLIC :: main_tracer_convec
  PUBLIC :: main_tracer_swap
#endif
  PUBLIC :: main_tracer_global_end    ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_free_memory   ! CALLED FROM messy_main_control_e5/c4
  !
  !PRIVATE :: setup_tracer_set_gp
  !PRIVATE :: setup_tracer_set_lg
  !PRIVATE :: setup_tracer_set_om
  !PRIVATE :: setup_tracer_set_ec
  !PRIVATE :: setup_tracer_set_cl
  !
  !PRIVATE :: main_tracer_read_nml_cpl

  ! SUBMODEL SMIL
  !
  ! CONVERSION ROUTINES (FAMILIES) TO BE USED IN SUBMODEL SMIL
  PUBLIC :: main_tracer_fconv_loc     ! CONVERT FAMILIES <-> TRACERS
  PUBLIC :: main_tracer_fconv_glb     ! CONVERT FAMILIES <-> TRACERS
  !
  !PRIVATE :: tracer_init ! initialise tracers
  !PRIVATE :: init_pvar   ! initialise prognostic variables of basemodel
#ifdef ICON
  !PRIVATE :: map_tracer_meta_i2m
  !PRIVATE :: map_tracer_meta_m2i
  !PRIVATE :: map_tracer_memory_icon
  !PRIVATE :: icon_tracer_integration
#endif

CONTAINS

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_setup

    ! TRACER MODULE ROUTINE (BASEMODEL INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, July 2002

    ! NOTES:
    ! - n_trcr_block must be read before initialize for COSMO

    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_tools,            ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_setup'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i

    ! select basemodel dependent rank of tracer index in internal memory
    TRRANK = _RI_TRRANK

    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL main_tracer_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi( &
            'main_tracer_read_nml_ctrl reported an error', substr)
    END IF
    CALL p_bcast(l_family, p_io)
    CALL p_bcast(l_pdef, p_io)
    DO i=1, NMAXTRACPROP
       CALL p_bcast(TPROP(i)%trset, p_io)
       CALL p_bcast(TPROP(i)%trlist, p_io)
       CALL p_bcast(TPROP(i)%caskname, p_io)
       CALL p_bcast(TPROP(i)%cont, p_io)
    END DO

    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL main_tracer_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi( &
            'main_tracer_read_nml_cpl reported an error', substr)
    END IF
    CALL p_bcast(l_tracer_init, p_io)
    CALL p_bcast(l_tracer_initfromrestart, p_io)
    CALL p_bcast(l_conv_lg2gp, p_io)
    CALL p_bcast(i_conv_lg2gp_mode, p_io)
    CALL p_bcast(r_conv_lg2gp_fill, p_io)
    CALL p_bcast(l_conv_lg2gp_mc, p_io)
    CALL p_bcast(n_trcr_block, p_io)
    !
    CALL p_bcast(l_conv_cl2gp, p_io)
    CALL p_bcast(i_conv_cl2gp_mode, p_io)
    CALL p_bcast(r_conv_cl2gp_fill, p_io)

    DO i=1, N_INI_MAX
       CALL p_bcast(ini_pvar(i)%basemodel, p_io)
       CALL p_bcast(ini_pvar(i)%pvar,      p_io)
       CALL p_bcast(ini_pvar(i)%inp,       p_io)
       CALL p_bcast(ini_pvar(i)%switch,    p_io)
    END DO

  END SUBROUTINE main_tracer_setup
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_initialize

    ! keep PDEF and FAMIILY CALLs in initialize because of timer_event_init

    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_initialize
    USE messy_main_tracer_pdef_bi,   ONLY: main_tracer_pdef_initialize

    IMPLICIT NONE

    IF (l_family) CALL main_tracer_family_initialize

    IF (l_pdef) CALL main_tracer_pdef_initialize

  END SUBROUTINE main_tracer_initialize
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_new_tracer(flag)

    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_new_tracer
    USE messy_main_channel_bi,       ONLY: REPR_UNDEF, GP_3D_MID &
                                         , GP_3D_MPIOM
!!#D clams +
    USE messy_main_channel_bi,       ONLY: REPR_LG_CLAMS
!!#D clams -
#if defined(ICON)
    USE messy_main_blather_bi,       ONLY: info_bi
    USE messy_main_tools,            ONLY: int2str
    USE messy_main_channel_mem,      ONLY: dom_current, n_dom
    USE messy_main_bmluse_bi,        ONLY: ntracer
#endif

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_new_tracer'
    INTEGER                     :: status
#if defined(ICON)
    INTEGER                     :: jg
    INTEGER                     :: zid
    CHARACTER(LEN=2)            :: pstr = '  '
    CHARACTER(LEN=3)            :: nstr = '   '
#endif

    ! INIT SWITCH

    L_GP = (GP_3D_MID  /= REPR_UNDEF)
!!#D clams +
    L_CL = (REPR_LG_CLAMS /= REPR_UNDEF)
!!#D clams -

    L_LG   = (NGCELL > 0)

    L_OM   = (GP_3D_MPIOM /= REPR_UNDEF)

    SELECT CASE(flag)
    CASE(1)
       !
       CALL start_message_bi(modstr,'SETUP TRACER SETS',substr)
       !
#if ! defined(ICON)
       CALL new_tracer_set(status, GPTRSTR, L_GP)
       CALL tracer_halt(substr, status)
#else
       ! SAVE ORIGINAL NUMBER OF ICON TRACERS, BEFORE MESSy TRACERS ARE ADDED
       ALLOCATE(l_ntrac_icon(n_dom))
       l_ntrac_icon(:) = ntracer

       ALLOCATE(L_GPTRSTR(n_dom))
       L_GPTRSTR(:) = ''
       DO jg=1, n_dom
          CALL int2str(pstr, jg)
          ! TRACER SET NAME
          L_GPTRSTR(jg)   = 'gp_D'//pstr
          CALL new_tracer_set(status, L_GPTRSTR(jg), .TRUE.)
          CALL tracer_halt(substr, status)
       END DO
#endif
       !
       ! added for clams tracer development
       CALL new_tracer_set(status, CLTRSTR, L_CL)
       CALL tracer_halt(substr, status)
       !
       CALL new_tracer_set(status, LGTRSTR, L_LG)
       CALL tracer_halt(substr, status)
       !
       CALL new_tracer_set(status, OMTRSTR, L_OM)
       CALL tracer_halt(substr, status)
       !
       CALL end_message_bi(modstr,'SETUP TRACER SETS',substr)
       !
#if defined(ICON)
       ! ... here, we are OUTSIDE the loop over patches!
       !
       CALL start_message_bi(modstr &
            ,'PLACE HOLDERS FOR ICON TRACERS',substr)
       !
       !     (for all patches)
       CALL map_tracer_meta_i2m(1)
       !
       CALL end_message_bi(modstr &
            ,'PLACE HOLDERS FOR ICON TRACERS',substr)
       !
#endif
       !
    CASE (2)
       !
       IF (l_family) CALL main_tracer_family_new_tracer
       !
#if defined(ICON)
       !
       ! UPDATE NUMBER OF TRACERS IN ICON
       ! Note: it seems that ICON can only manage the same number of tracers
       !       for all patches, at least ntracer does not depend on the patch
       CALL start_message_bi(modstr, 'UPDATE NUMBER OF TRACERS IN ICON' &
            , substr)
       !
       DO jg=1, n_dom
          CALL get_tracer_set_id(status, L_GPTRSTR(jg), zid)
          CALL tracer_halt(substr//'(2)', status)
          !
          !ntracer = trset(zid)%ntrac
          ntracer = MAX(trset(zid)%ntrac, ntracer)
          !
          CALL int2str(pstr, jg)
          CALL int2str(nstr, ntracer,' ')
          CALL info_bi(' new number of tracers: '//nstr//' (patch='//pstr//')' &
               , substr)
       END DO
       !
       CALL end_message_bi(modstr, 'UPDATE NUMBER OF TRACERS IN ICON' &
            , substr)
#endif
       !
    CASE(3)
       !
#ifndef ICON
       CALL start_message_bi(modstr,'OVERWRITE TRACER PROPERTIES',substr)
       CALL set_tracer_properties(status, lprint=p_parallel_io)
       CALL tracer_halt(substr, status)
       CALL end_message_bi(modstr,'OVERWRITE TRACER PROPERTIES',substr)
       !
       CALL start_message_bi(modstr,'SHOW TRACER SETS',substr)
       !
       IF (p_parallel_io) CALL print_tracer_set
       !
       CALL end_message_bi(modstr,'SHOW TRACER SETS',substr)
#endif
       !
    CASE DEFAULT
       !
       CALL error_bi( 'UNKNOWN FLAG !', substr)
       !
    END SELECT

  END SUBROUTINE main_tracer_new_tracer
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_init_memory(flag)

    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_tools,            ONLY: remap_bounds
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
#ifdef ECHAM5
    USE messy_main_channel_bi,       ONLY: LG_ATTILA, GP_3D_MPIOM
#endif
!!#D clams +
    USE messy_main_channel_bi,       ONLY: REPR_LG_CLAMS
!!#D clams -
    USE messy_main_channel_tracer,   ONLY: create_tracer_channels
    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_init_mem
    USE messy_main_tracer_pdef_bi,   ONLY: main_tracer_pdef_init_mem
#ifdef COSMO
    USE messy_main_data_bi,          ONLY: L_IS_CHILD, l2tls
#endif
#ifdef ICON
    USE messy_main_tools,         ONLY: int2str
    USE messy_main_channel_mem,   ONLY: dom_current
#endif

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)           :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_init_memory'
    INTEGER :: status
#ifdef ICON
    CHARACTER(LEN=2) :: pstr = '  '
#endif

    SELECT CASE(flag)
    CASE (1) !==============================================================
       !
       CALL start_message_bi(modstr,'SETUP TRACER MEMORY',substr)

       CALL setup_tracer_set_gp

       ! INITALIZE POINTER (messy_main_tracer_mem_bi.90)
       IF (L_GP) THEN
          xt    => pxt(:,:,:,:,1)
          xtm1  => pxtm1(:,:,:,:,1)
          xtte  => pxtte(:,:,:,:,1)
#if ! (defined(COSMO) || defined(ICON))
          xtf   => pxtf(:,:,:,:,1)
#else
#ifndef ICON
          IF (l2tls) THEN
             xtf   => pxtf(:,:,:,:,1)
          END IF
#endif
#endif
       ENDIF

#ifdef COSMO
       ! correct pointer from above and allocate pointer array xt_array
       IF (l2tls) THEN
          NULLIFY(xtf)
          xt_bd => pxtf(:,:,:,:,1:2)
       ELSE
          xtf   => pxtf(:,:,:,:,1)
          xt_bd => pxtf(:,:,:,:,2:3)
       ENDIF
       ! the boundary data is defined as place 1 and 2 of
       ! extended tracer memory
       xt_bd(:,:,:,:,:) = -999._dp
#endif

#if defined(ECHAM5)
       IF ( l_conv_cl2gp .AND. (.NOT. L_CL ) ) THEN
          IF (p_parallel_io) THEN
             WRITE(*,*) ' L_CONV_CL2GP = T in namelist will be ignored,'
             WRITE(*,*) ' since CLAMS is NOT running!'
          END IF
          l_conv_cl2gp = .FALSE.
       END IF
       CALL setup_tracer_set_cl

       ! initalize CLaMS TRACER pointer (messy_main_tracer_mem_bi.90)
       if (l_cl) then
          xt_c    => pxt_c(:,:,:,:,1)
          xtm1_c  => pxt_c(:,:,:,:,1)  !!! same as xt for CLamS !!!
          xtte_c  => pxtte_c(:,:,:,:,1)
          IF (l_conv_cl2gp) xt_clgp => pxt_clgp(:,:,:,:,1)
       endif

       ! ------------

       IF ( l_conv_lg2gp .AND. (.NOT. L_LG ) ) THEN
          IF (p_parallel_io) THEN
             WRITE(*,*) ' L_CONV_LG2GP = T in namelist will be ignored,'
             WRITE(*,*) ' since no LAGRANGIAN scheme is running!'
          END IF
          l_conv_lg2gp = .FALSE.
       END IF

       CALL setup_tracer_set_lg

       ! INITALIZE POINTER (messy_main_tracer_mem_bi.90)
       IF (L_LG) THEN
          xt_a => pxt_a(:,:,:,:,1)
          xtm1_a => pxtm1_a(:,:,:,:,1)
          xtte_a => pxtte_a(:,:,:,:,1)
          xtf_a  => pxtf_a(:,:,:,:,1)
          !
          IF (l_conv_lg2gp) xt_lggp => pxt_lggp(:,:,:,:,1)
       END IF
       !

       CALL setup_tracer_set_om
       ! INITALIZE POINTER (messy_main_tracer_mem_bi.90)
       IF (L_OM) THEN
          xt_om => pxt_om(:,:,:,:,1)
       END IF
#endif
       !
       CALL end_message_bi(modstr,'SETUP TRACER MEMORY',substr)
       !
    CASE (2) !==============================================================
       !
       CALL start_message_bi(modstr,'TRACER MEMORY CHANNEL COUPLING',substr)
       !
       ! COUPLE TRACER MEMORY TO CHANNELS
#ifndef COSMO
       CALL create_tracer_channels(status, GPTRSTR, gp_channel, GP_3D_MID)
       CALL channel_halt(substr, status)
#else
       IF (L_IS_CHILD) THEN
          CALL create_tracer_channels(status, GPTRSTR, gp_channel, GP_3D_MID, 2)
          CALL channel_halt(substr, status)
       ELSE
          CALL create_tracer_channels(status, GPTRSTR, gp_channel, GP_3D_MID)
          CALL channel_halt(substr, status)
       ENDIF
#endif
       !
#if defined(ECHAM5)
!!#D clams +
       ! creating CLaMS tracer channels
       CALL create_tracer_channels(status, &
           CLTRSTR,   cl_channel  , REPR_LG_CLAMS)
       CALL channel_halt(substr, status)
       IF (l_conv_cl2gp) THEN
          CALL create_tracer_channels(status, &
               CLGPTRSTR, clgp_channel, GP_3D_MID)
          CALL channel_halt(substr, status)
       END IF
!!#D clams -

       CALL create_tracer_channels(status, LGTRSTR, lg_channel, LG_ATTILA)
       CALL channel_halt(substr, status)
       !
       IF (l_conv_lg2gp) THEN
          CALL create_tracer_channels(status, LGGPTRSTR, &
               lggp_channel, GP_3D_MID)
          CALL channel_halt(substr, status)
       END IF
       !
       CALL create_tracer_channels(status, OMTRSTR, om_channel, GP_3D_MPIOM)
       CALL channel_halt(substr, status)
#endif
       !
#if ! defined(ICON)
       IF (l_pdef) CALL main_tracer_pdef_init_mem
       !
       ! setting meta information of family-members to fraction
       ! (for advection initialization)
       IF (l_family) CALL main_tracer_family_init_mem
       !
#endif
!!$#if defined(ICON)
!!$       !
!!$       CALL map_tracer_memory_icon(2)
!!$       !
!!$#endif
       !
       CALL end_message_bi(modstr,'TRACER MEMORY CHANNEL COUPLING',substr)
       !
#if defined(ICON)
       !
    CASE(11) ! REDIRECT ICON TRACER MEMORY ------------------------------
       !
       ! ... here, we are INSIDE a loop over patches
       !
       CALL int2str(pstr, dom_current)
       !
       CALL start_message_bi(modstr &
            ,'MAP ICON TRACER MEMORY TO MESSy TRACER MEMORY (patch='//&
            &pstr//')',substr)
       ! overwrite ntracer in ICON, deallocate ICON tracer pointers and
       ! point them to the corresponding MESSy TRACER memory
       ! (only current patch)
!qqq
       CALL map_tracer_memory_icon(1)
       !
       CALL end_message_bi(modstr &
            ,'MAP ICON TRACER MEMORY TO MESSy TRACER MEMORY (patch='//&
            &pstr//')',substr)
       !
    CASE(12) ! UPDATE MESSY TRACER META INFORMATION OF ICON TRACERS -----
       !
       ! ... here, we are INSIDE a loop over patches
       !
       CALL int2str(pstr, dom_current)
       !
       CALL start_message_bi(modstr &
            ,'UPDATE MESSY TRACER META INFORMATION OF ICON TRACERS (patch='//&
            &pstr//')',substr)
       !
       CALL map_tracer_meta_i2m(2)
       !
       CALL end_message_bi(modstr &
            ,'UPDATE MESSY TRACER META INFORMATION OF ICON TRACERS (patch='//&
            &pstr//')',substr)
       !
       ! - OVERWRITE TRACER PROPERTIES VIA NAMELIST ----------------
       ! NOTE: This is now done multiple times for each patch,
       !       since it is called from within a loop over patches ...
       CALL start_message_bi(modstr,'OVERWRITE TRACER PROPERTIES',substr)
       !
       CALL set_tracer_properties(status, lprint=p_parallel_io &
            , dom_id = dom_current)
       CALL tracer_halt(substr, status)
       !
       CALL end_message_bi(modstr,'OVERWRITE TRACER PROPERTIES',substr)
       !
       ! - MAP MESSy TRACER META INFORMATION TO ICON TRACER META INFORMATION
       !
       CALL start_message_bi(modstr &
            ,'MESSy TRACER TO ICON MAPPING (patch='//pstr//')',substr)
       !
       ! map MESSy TRACER meta information to ICON tracer meta information
       ! (only current patch)
       CALL map_tracer_meta_m2i
       !
       CALL end_message_bi(modstr &
            ,'MESSy TRACER TO ICON MAPPING (patch='//pstr//')',substr)
       !
    CASE(3) !==============================================================
       !
       IF (l_pdef) CALL main_tracer_pdef_init_mem
       !
       ! setting meta information of family-members to fraction
       ! (for advection initialization)
       IF (l_family) CALL main_tracer_family_init_mem
       !
       CALL start_message_bi(modstr,'SHOW TRACER SETS',substr)
       !
       IF (p_parallel_io) CALL print_tracer_set
       !
       CALL end_message_bi(modstr,'SHOW TRACER SETS',substr)
       !
#endif
       !
    CASE DEFAULT
       !
       CALL error_bi('UNKNOWN FLAG !', substr)
       !
    END SELECT

  END SUBROUTINE main_tracer_init_memory
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_init_coupling(flag)

    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_init_cpl
    USE messy_main_tracer_pdef_bi,   ONLY: main_tracer_pdef_init_cpl
    USE messy_main_blather_bi,       ONLY: error_bi

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: flag

    CHARACTER(LEN=*), PARAMETER :: substr='main_tracer_init_coupling'

    SELECT CASE(flag)

    CASE (1)
       IF (l_family) CALL main_tracer_family_init_cpl

       IF (l_pdef) CALL main_tracer_pdef_init_cpl
#ifdef ICON
    CASE (2)
       CALL map_tracer_memory_icon(2)
#endif
    CASE DEFAULT
       CALL error_bi('UNKNOWN FLAG ! ',substr)
    END SELECT

  END SUBROUTINE main_tracer_init_coupling
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_init_tracer(flag)

    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: error_bi, info_bi
    USE messy_main_data_bi,          ONLY: eps
    USE messy_main_timer,            ONLY: lresume
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_tracer,           ONLY: NSETID, STRLEN_TRSET &
                                         , get_tracer_set       &
                                         , I_MMD_INIT, ON
    USE messy_main_channel,          ONLY: get_channel_object_info
#ifdef ICON
    USE messy_main_tools,            ONLY: split_name_domain
#endif

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, NULL, TRIM

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'main_tracer_init_tracer'
    INTEGER                      :: status
    TYPE(t_trinfo_list), POINTER :: til => NULL()
    CHARACTER(LEN=STRLEN_TRSET)  :: setname = ''
    REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: zpxt   => NULL()
    REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: zpxtm1 => NULL()
    LOGICAL :: linit, linit_m1
    LOGICAL :: l_init
    LOGICAL :: linit_mmd
    INTEGER :: i, ntrac
    LOGICAL :: l_enabled
#ifdef ICON
    INTEGER                     :: domain
#endif
    CHARACTER(LEN=STRLEN_TRSET) :: sname = ''

    SELECT CASE(flag)
    CASE(1)
       !
       CALL start_message_bi(modstr,'CHECK TRACER INIT FROM RESTART',substr)
       !
    CASE(4)
       !
       CALL start_message_bi(modstr,'INITIALISE TRACER VIA TRACER_INIT',substr)
       IF (l_tracer_init) THEN
          CALL tracer_init
       ELSE
          CALL info_bi(' ... skipped (l_tracer_init=F in CPL)', substr)
       END IF
       CALL end_message_bi(modstr,'INITIALISE TRACER VIA TRACER_INIT',substr)
       RETURN
       !
    CASE(2)
       CALL start_message_bi(modstr,'CHECK TRACER INIT BY TRACER_INIT',substr)
       !
    CASE(3)
       !
       CALL start_message_bi(modstr,'DIAGNOSE TRACER INITIALIZATION',substr)
       !
    CASE DEFAULT
       !
       CALL error_bi('UNKNOWN FLAG !', substr)
       !
    END SELECT

    set_loop: DO i=1, NSETID
       !
       CALL get_tracer_set(status, i, setname=setname &
            , trlist=til, xt=zpxt, xtm1=zpxtm1 &
            , ntrac=ntrac, l_init=l_init       &
            , l_enable=l_enabled)
       CALL tracer_halt(substr, status)
       IF (p_parallel_io) THEN
          WRITE(*,*) '*** TRACER SET '//TRIM(setname)//' *** enabled = '&
               , l_enabled
       END IF

       IF (.NOT. l_init)    CYCLE
       IF (.NOT. l_enabled) CYCLE
       IF (ntrac <= 0)      CYCLE

       tracer_loop: DO
          IF (.NOT.ASSOCIATED(til)) EXIT

          SELECT CASE(flag)
          CASE(1)
             !
             ! CHECK IF TRACER HAS BEEN INITIALIZED VIA RESTART
             ! -> CHANNEL OBJECTS FOR X and X_m1 have restart_read = .true.
#ifdef ICON
             CALL split_name_domain(status, setname, sname, domain)
             sname = modstr//'_'//TRIM(sname)
#else
             sname = modstr//'_'//TRIM(setname)
#endif
             CALL get_channel_object_info(status, TRIM(sname)           &
                  , TRIM(til%info%ident%fullname), lrestart_read=linit  &
#ifdef ICON
                  , dom_id = domain &
#endif
                  )
             CALL channel_halt(substr, status)

             IF (ASSOCIATED(zpxtm1)) THEN
                CALL get_channel_object_info(status, TRIM(sname)//'_m1' &
                     , TRIM(til%info%ident%fullname), lrestart_read=linit_m1 &
#ifdef ICON
                     , dom_id = domain &
#endif
                     )
                CALL channel_halt(substr, status)
             ELSE
                linit_m1 = .TRUE.
             END IF

             ! CHECK IF gridpoint tracer was initialised by MMD
             linit_mmd = til%info%meta%cask_i(I_MMD_INIT) == ON
             ! ... or by BASEMODEL
             til%info%meta%linit = til%info%meta%linit .OR. &
                  (linit .AND. linit_m1) .OR. linit_mmd

             IF (til%info%meta%linit) THEN
                IF (p_parallel_io) THEN
#ifndef I2CINC
                   WRITE(*,*) 'TRACER '//TRIM(til%info%ident%fullname)//' (', &
                        til%info%ident%idx,                  &
                        ') WAS INITIALIZED FROM RESTART-FILE!'
#else
                   IF (linit_mmd .AND. (.NOT. linit) ) THEN
                   WRITE(*,*) 'TRACER '//TRIM(til%info%ident%fullname)//' (', &
                        til%info%ident%idx,                  &
                        ') WAS INITIALIZED BY MMD!'
                   ELSE
                   WRITE(*,*) 'TRACER '//TRIM(til%info%ident%fullname)//' (', &
                        til%info%ident%idx,                  &
                        ') WAS INITIALIZED FROM RESTART-FILE!'
                   ENDIF
#endif
                END IF
                ! CHECK FOR lforce_init (FORCE TRACER_INIT OR vini)
                IF (til%info%meta%cask_i(I_force_init)==ON) THEN
                   IF (p_parallel_io) THEN
                      WRITE(*,*) ' force_init -> ', &
                           'RESTART INITIALIZATION WILL BE IGNORED'
                   END IF
                   til%info%meta%linit = .false.
                END IF
             ELSE
                IF (p_parallel_io) THEN
                   WRITE(*,*) 'TRACER '//TRIM(til%info%ident%fullname)//' (', &
                        til%info%ident%idx,                  &
                        ') WAS NOT INITIALIZED FROM RESTART-FILE ...'
                END IF
             END IF
             !
          CASE(2)
             !
             ! CHECK IF TRACER HAS BEEN INITIALIZED IN SOME WAY
             IF (til%info%meta%linit) THEN
                IF (p_parallel_io) THEN
                   WRITE(*,*) 'TRACER '//TRIM(til%info%ident%fullname)//' (', &
                        til%info%ident%idx,                  &
                        ') ALREADY INITIALIZED ... !'
                END IF
             ELSE
                IF (p_parallel_io) THEN
                   WRITE(*,*) 'TRACER '//TRIM(til%info%ident%fullname)//' (', &
                        til%info%ident%idx,                  &
                        ') WILL BE SET TO ', til%info%meta%cask_R(R_vini)
                END IF

                zpxt( _RI_XYZN_(:,:,:,til%info%ident%idx),:) = &
                     til%info%meta%cask_R(R_vini)

                IF (ASSOCIATED(zpxtm1)) THEN
                   IF (lresume) &
                        zpxtm1(_RI_XYZN_(:,:,:,til%info%ident%idx),:) = &
                        (1._DP -eps) &
                        * zpxt(_RI_XYZN_(:,:,:,til%info%ident%idx),:)
                END IF

                til%info%meta%linit = .TRUE.
             END IF
             !
          CASE (3)
             !
             ! NOTHING TO DO
          END SELECT

          til => til%next
       END DO tracer_loop
       !
    END DO set_loop


    SELECT CASE(flag)
    CASE(1)
       !
       CALL end_message_bi(modstr,'CHECK TRACER INIT FROM RESTART',substr)
       !
    CASE(2)
       !
       CALL end_message_bi(modstr,'CHECK TRACER INIT BY TRACER_INIT',substr)
       !
    CASE(3)
       !
       ! --- DIAGNOSTIC OUTPUT
       !
       IF (p_parallel_io) CALL print_tracer_set_val
       !
       CALL end_message_bi(modstr,'DIAGNOSE TRACER INITIALIZATION',substr)
       !
    CASE DEFAULT
       !
       CALL error_bi('UNKNOWN FLAG !', substr)
       !
    END SELECT

  END SUBROUTINE main_tracer_init_tracer
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_global_start

    USE messy_main_tracer_pdef_bi,     ONLY: main_tracer_pdef_global_start

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_global_start'

#if defined(ECHAM5)
    ! LAGRANGIAN TRACERS ATTILA
    IF (L_LG) THEN
       qxt_a    => xt_a(_RI_XYZN_(:,1,1,:))
       qxtte_a  => xtte_a(_RI_XYZN_(:,1,1,:))
       qxtm1_a  => xtm1_a(_RI_XYZN_(:,1,1,:))
       qxtf_a   => xtf_a(_RI_XYZN_(:,1,1,:))
    ENDIF

    ! LAGRANGIAN TRACERS CLaMS
    IF (L_CL) THEN
       qxt_c    => xt_c(_RI_XYZN_(:,1,1,:))
       qxtte_c  => xtte_c(_RI_XYZN_(:,1,1,:))
       qxtm1_c  => xtm1_c(_RI_XYZN_(:,1,1,:))
    ENDIF
#endif
       !
       IF (l_pdef) CALL main_tracer_pdef_global_start

#if defined(ICON) && !defined(MESSYTENDENCY)
       ! reset MESSy tracer tendencies for current patch
       xtte(:,:,:,:) = 0.0_dp
#endif

  END SUBROUTINE main_tracer_global_start
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_beforeadv

    USE messy_main_tracer_family_bi,   ONLY: main_tracer_family_beforeadv

    IMPLICIT NONE

       ! TYPE-1: t2f
       ! TYPE-2: summation (GPTRSTR)
       IF (l_family) CALL main_tracer_family_beforeadv

  END SUBROUTINE main_tracer_beforeadv
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_afteradv

    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_afteradv

    IMPLICIT NONE

    IF (l_family) CALL main_tracer_family_afteradv

  END SUBROUTINE main_tracer_afteradv
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_local_start

    USE messy_main_grid_def_mem_bi,  ONLY: jrow
#ifdef COSMO
    USE messy_main_data_bi,          ONLY: l2tls
#endif

    IMPLICIT NONE

    ! SET POINTERS FOR WITHIN LOCAL LOOP

    ! GRIDPOINT TRACERS
    IF (L_GP) THEN
       qxt    => xt(_RI_XYZN_(:,jrow,:,:))
       qxtte  => xtte(_RI_XYZN_(:,jrow,:,:))
       qxtm1  => xtm1(_RI_XYZN_(:,jrow,:,:))
#if defined(ECHAM5) || defined(MBM_CLAMS) || defined(VERTICO) || defined(MBM_MPIOM)
       qxtf   => xtf(_RI_XYZN_(:,jrow,:,:))
#endif
# if defined(COSMO)
       IF (.NOT. l2tls) THEN
          qxtf   => xtf(_RI_XYZN_(:,jrow,:,:))
       END IF
#endif
    END IF

#if defined(ECHAM5)
    IF (L_LG) THEN
       IF (l_conv_lg2gp) qxt_lggp => xt_lggp(_RI_XYZN_(:,jrow,:,:))
    END IF

    IF (L_CL) THEN
       IF (l_conv_cl2gp) qxt_clgp => xt_clgp(_RI_XYZN_(:,jrow,:,:))
    END IF
#endif

  END SUBROUTINE main_tracer_local_start
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_global_end

#if ! ( defined(BLANK) || defined(MBM_CLAMS) || defined(VERTICO) || defined(MBM_MPIOM) )
    USE messy_main_grid_def_bi,      ONLY: grmass=>grmassdry
#else
    USE messy_main_blather_bi,       ONLY: warning_bi
#endif
    USE messy_main_tracer_pdef_bi,   ONLY: main_tracer_pdef_global_end
    USE messy_main_tracer_pdef,      ONLY: tracpdef_airmass
#if defined(ECHAM5)
!!#D attila +
    USE messy_attila_mem,            ONLY: AMCELL
!!#D attila -
!!#D mpiom +
    USE messy_mpiom_e5,              ONLY: omvol
!!#D mpiom -
#endif

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_global_end'

#ifdef ICON
    ! add all MESSy tracer tendencies to actual tracers
    CALL icon_tracer_integration
#endif

    IF (l_pdef) THEN
#if ! ( defined(BLANK) || defined(MBM_CLAMS) || defined(VERTICO) || defined(MBM_MPIOM) )
       !
       CALL tracpdef_airmass(GPTRSTR, grmass)
       !
#if defined(ECHAM5)
!!#D attila +
       IF (L_LG) CALL tracpdef_airmass(LGTRSTR, AMCELL)
       IF (l_conv_lg2gp) CALL tracpdef_airmass(LGGPTRSTR, grmass)
!!#D attila -
       !
!!#D mpiom +
       IF (L_OM) CALL tracpdef_airmass(OMTRSTR, omvol)
!!#D mpiom -
#endif
       !
       CALL main_tracer_pdef_global_end
       !
#else
#ifdef BLANK
       CALL warning_bi(substr, 'pdef for BLANK not possible')
#endif
#ifdef MBM_MPIOM
       CALL warning_bi(substr, 'pdef for MBM_MPIOM not possible')
#endif
#ifdef MBM_CLAMS
       CALL warning_bi(substr, 'pdef for CLAMS not yet possible')
#endif
#ifdef VERTICO
       CALL warning_bi(substr, 'pdef for VERTICO not possible')
#endif
       CALL warning_bi(substr, 'define grmass first')
#endif
    END IF

  END SUBROUTINE main_tracer_global_end
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
#ifdef ICON
  SUBROUTINE main_tracer_swap(jg)

    USE messy_main_bmluse_bi,  ONLY: p_nh_state, nnow_rcf, nnew_rcf

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: jg

    ! LOCAL
    INTEGER             :: jt

    ! copy new integrated value as start value for next time step
    pxtm1(:,:,:,:,1) = pxt(:,:,:,:,1)
    !qqq do we need nnew/ nnow or nnew_rcf /nnow_rcf here?
    p_nh_state(jg)%prog(nnow_rcf(jg))%tracer => pxtm1(:,:,:,:,1)
    p_nh_state(jg)%prog(nnew_rcf(jg))%tracer => pxt(:,:,:,:,1)

    DO jt = 1, SIZE(p_nh_state(jg)%prog(nnew_rcf(jg))%tracer_ptr)
       p_nh_state(jg)%prog(nnow_rcf(jg))%tracer_ptr(jt)%p_3d => pxtm1(_RI_XYZN_(:,:,:,jt),1)
       p_nh_state(jg)%prog(nnew_rcf(jg))%tracer_ptr(jt)%p_3d => pxt(_RI_XYZN_(:,:,:,jt),1)
    END DO

  END SUBROUTINE main_tracer_swap
#endif
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_free_memory

    USE messy_main_tracer_pdef_bi,   ONLY: main_tracer_pdef_free_mem
    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_free_mem
#if defined(ICON)
    USE messy_main_bmluse_bi,        ONLY: prm_nwp_tend
#endif

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_free_memory'
    INTEGER :: status
#ifdef ICON
    INTEGER :: jg, jb, ic
#endif

    CALL start_message_bi(modstr,'FREE TRACER MEMORY',substr)

#if ! defined(ICON)
    CALL clean_tracer_set(status, GPTRSTR)
    CALL tracer_halt(substr, status)
#else
    DO jg=1, SIZE(L_GPTRSTR)

       CALL clean_tracer_set(status, L_GPTRSTR(jg))
       CALL tracer_halt(substr, status)

       ! free memory for convective tracer tendencies
       ! Note: this must be called before destruct_nwp_phy_state, which
       !       is called from mo_atmo_nonhydrostatic!
       IF (ALLOCATED(prm_nwp_tend)) THEN
          DO ic=1, SIZE(prm_nwp_tend(jg)%conv_tracer_tend,2)
             DO jb=1, SIZE(prm_nwp_tend(jg)%conv_tracer_tend,1)
!                DEALLOCATE(prm_nwp_tend(jg)%conv_tracer_tend(jb, ic)%ptr)
                NULLIFY(prm_nwp_tend(jg)%conv_tracer_tend(jb, ic)%ptr)
             END DO
          END DO
          DO ic=1, SIZE(prm_nwp_tend(jg)%turb_tracer_tend,2)
             DO jb=1, SIZE(prm_nwp_tend(jg)%turb_tracer_tend,1)
!                DEALLOCATE(prm_nwp_tend(jg)%conv_tracer_tend(jb, ic)%ptr)
                NULLIFY(prm_nwp_tend(jg)%turb_tracer_tend(jb, ic)%ptr)
             END DO
          END DO
       END IF

    END DO
#endif

#if defined(ECHAM5)
    ! CLaMS tracer set added, cleaning here
    CALL clean_tracer_set(status, CLTRSTR)
    CALL tracer_halt(substr, status)
    IF (l_conv_cl2gp) THEN
       CALL clean_tracer_set(status, CLGPTRSTR)
       CALL tracer_halt(substr, status)
    END IF
#endif

#if defined(ECHAM5)
    CALL clean_tracer_set(status, LGTRSTR)
    CALL tracer_halt(substr, status)
    !
    IF (l_conv_lg2gp) THEN
       CALL clean_tracer_set(status, LGGPTRSTR)
       CALL tracer_halt(substr, status)
    END IF
#endif

    IF (l_pdef) CALL main_tracer_pdef_free_mem

    IF (l_family) CALL main_tracer_family_free_mem

#if defined(ICON)
       DEALLOCATE(L_GPTRSTR) ; NULLIFY(L_GPTRSTR)
       !DEALLOCATE(l_gp_channel) ; NULLIFY(l_gp_channel)
       DEALLOCATE(l_ntrac_icon) ; NULLIFY(l_ntrac_icon)
       DO jg=1, SIZE(grmass_conv)
          DEALLOCATE(grmass_conv(jg)%ptr);    NULLIFY(grmass_conv(jg)%ptr)
          DEALLOCATE(grmassdry_conv(jg)%ptr); NULLIFY(grmassdry_conv(jg)%ptr)
       END DO
       DEALLOCATE(grmass_conv);    NULLIFY(grmass_conv)
       DEALLOCATE(grmassdry_conv); NULLIFY(grmassdry_conv)
#endif

    CALL end_message_bi(modstr,'FREE TRACER MEMORY',substr)

  END SUBROUTINE main_tracer_free_memory
  ! -------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE TRACER_INIT

#if defined(ECHAM5) || defined(COSMO) || defined(CESM1)

    ! INITIALIZES TRACER FIELDS FROM netCDF FILES
    ! (USING IMPORT_GRID) AS DEFINED IN NAMELIST FILE tracer.nml'
    !
    ! Author: Patrick Joeckel, MPICH, Mainz, November 2002
    !         Patrick Joeckel, MPICH, Mainz, March    2004
    !         Patrick Joeckel, DLR-IPA, Oberpfaffenhofen, December 2009
    !         Astrid  Kerkweg, UNI-Mz, Mainz, August 2013

    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_io, p_bcast  &
                                         , dcg, scatter_gp               &
                                         , p_pe, p_nprocs, p_all_comm    &
                                         , p_send, p_recv
    USE messy_main_blather_bi,      ONLY: start_message_bi, end_message_bi &
                                        , error_bi, warning_bi
    USE messy_main_data_bi,          ONLY: eps
#if defined(COSMO) || defined(MESSYDWARF)
    USE messy_main_data_bi,          ONLY: basemodstr => modstr
#endif
    USE messy_main_grid_def_mem_bi,  ONLY: nlev, nproma, ngpblks        &
#ifndef COSMO
                                         , nlon, ngl                    &
#endif
                                         , BASEGRID_ID
    USE messy_main_import_grid,      ONLY: RG_CTRL, RG_NML              &
                                         , NML_NEXT   &
                                         , RG_SCAN, RG_STATUS           &
                                         , RGSTAT_STOP, READ_CONTROL
    USE messy_main_import_grid_par,  ONLY: RGSTAT_NULL, INIT_PIMPGRID   &
                                         , INIT_PARALLEL, pe_list       &
                                         , nproc_work
    USE messy_main_grid,             ONLY: t_geohybgrid, INIT_GEOHYBGRID &
                                         , grid_error
    USE messy_main_grid,             ONLY: LOCATE_GEOHYBGRID           &
                                         , INIT_GEOHYBGRID
    USE messy_main_grid_netcdf,      ONLY: RGMSG, ERRMSG, RGMLVL, RGMLW     &
                                         , RGMLWC, INIT_NCVAR &
                                         , t_ncvar
    USE messy_main_grid_tools,       ONLY: RGTOOL_CONVERT
    USE messy_main_grid_trafo_nrgd,  ONLY: REGRID_CONTROL
    USE messy_main_channel_error_bi, ONLY: CHANNEL_HALT
    USE messy_main_channel_bi,       ONLY: main_channel_read_restart
    USE messy_main_channel_mem,      ONLY: dom_curid
    USE messy_main_channel_repr,     ONLY: repr_def_axes
    USE messy_main_channel,          ONLY: GET_CHANNEL_OBJECT_INFO
    USE messy_main_constants_mem,    ONLY: STRLEN_MEDIUM, STRLEN_ULONG
    USE messy_main_timer,            ONLY: lresume, lstart
    USE messy_main_tracer,           ONLY: NSETID, STRLEN_TRSET, get_tracer_set&
                                         , t_trinfo_list, TR_NEXIST
    USE messy_main_tools,            ONLY: find_next_free_unit

#if defined(ECHAM5)
    ! ... LGTRSTR
!!#D attila +
    USE messy_attila_tools_e5,    ONLY: gp2lg_e5
!!#D attila -
!!#D clams +
    ! ... CLTRST
    USE messy_clams_tools_e5,     ONLY: gp2cl_e5
!!#D clams -
#endif

#if defined(COSMO) || (defined(CESM1) && defined(HOMMESE))
    USE messy_main_grid_trafo,      ONLY: COMPLETE_GEOHYBGRID
    USE messy_main_grid_trafo_scrp, ONLY: CALC_SCRIPDATA, CALC_SCRIP_WEIGHTS
    USE messy_main_grid_trafo_scrp, ONLY: CONSTRUCT_VERTICAL_AXIS
    USE messy_main_grid_trafo_scrp, ONLY: SCRIP_CONTROL
    USE messy_main_grid_trafo_scrp, ONLY: INTERPOL_GEOHYBGRID_PS
    USE messy_main_grid_trafo_scrp, ONLY: t_scrip_data
    USE messy_main_grid_trafo_nrgd_base, ONLY: SNREGRID
    USE messy_main_grid,            ONLY: GEOHYBGRID_UPDATE_PRESS
    USE messy_main_grid,            ONLY: NEW_GEOHYBGRID
    USE messy_main_grid_netcdf,     ONLY: RGMLI
    USE messy_main_grid_netcdf,     ONLY: QDEF_NCVAR
    USE messy_main_grid_tools,      ONLY: RGTOOL_CONVERT, RGTOOL_CONVERT_DAT2VAR
    USE messy_main_grid_tools,      ONLY: RGTOOL_CONVERT_DAT2PREDEFVAR
    USE messy_main_channel,         ONLY: t_chaobj_cpl
    USE messy_main_channel,         ONLY: GET_CHANNEL_OBJECT
    USE messy_main_channel_repr,    ONLY: GET_REPRESENTATION_INFO
    USE messy_main_tools,           ONLY: lcase
#endif

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, NULL, SIZE, TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'tracer_init'
    INTEGER                      :: n       ! set counter
    INTEGER                      :: ntrac   ! number of tracers
    CHARACTER(LEN=STRLEN_TRSET)  :: setname ! tracer set name
    LOGICAL                      :: l_init  ! initialize this set ?
    TYPE(t_trinfo_list), POINTER :: til => NULL()
    TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti  => NULL()
    !
    INTEGER                      :: status  ! error status
    LOGICAL                      :: lskip   ! skip, if no uninitialized tracers
    INTEGER                      :: iunit   ! fortran unit for input file
    LOGICAL                      :: lex     ! file exists ?
    !
    ! FOR REGRIDDING TO GLOBAL FIELD (ON I/O PE)
    REAL(DP), DIMENSION(:,:,:), POINTER :: zin => NULL()
    !
    ! decomposed data field (local)
    REAL(DP), DIMENSION(:,:,:), POINTER :: zinl => NULL()
#ifdef ECHAM5
    ! Lagrange
    REAL(DP), DIMENSION(:),     POINTER :: zxt_lg => NULL()
    REAL(DP), DIMENSION(:),     POINTER :: zxt_cl => NULL()
#endif
    !
    ! GRID
    TYPE (t_ncvar), DIMENSION(:), POINTER :: var => NULL() ! list of variables
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: dat => NULL()
    CHARACTER(len=2*STRLEN_MEDIUM + 1)    :: cevar ! object of evar
    CHARACTER(len=STRLEN_MEDIUM)          :: basename
    CHARACTER(len=STRLEN_MEDIUM)          :: subname
    INTEGER                               :: isizev ! SIZE(evar)
    INTEGER                               :: i      ! species counter
    INTEGER, DIMENSION(:), ALLOCATABLE    :: nu, nr, ni ! counter
    INTEGER                               :: jt         ! tracer index in set
    INTEGER                               :: root_pe
    INTEGER                               :: work_loop
    LOGICAL,  SAVE                        :: lint    ! use input time ?
    INTEGER, DIMENSION(:), POINTER        :: RGT => NULL() ! regridding type
    TYPE (t_ncvar), DIMENSION(:), POINTER :: ovar  => NULL() ! list of variables
    TYPE (t_ncvar), DIMENSION(:), POINTER :: sovar => NULL() ! list of variables
    TYPE (t_geohybgrid)                   :: igrid   ! output grid info
    TYPE (t_geohybgrid)                   :: ogrid   ! output grid info
    INTEGER                               :: tc, tmin
    INTEGER,           PARAMETER          :: t_undef = -99
    LOGICAL                               :: ldompar = .FALSE.
    LOGICAL                               :: lvarpar = .FALSE.
    LOGICAL                               :: l_proc  = .FALSE.
    LOGICAL                               :: lpres   = .FALSE.
    LOGICAL                               :: lwork   = .FALSE.
    INTEGER                               :: ix
    INTEGER                     :: np
    CHARACTER(LEN=STRLEN_ULONG) :: infostr = ''
    LOGICAL :: linit
#if defined(COSMO) || (defined(CESM1) && defined(HOMMESE))
    CHARACTER(LEN=4)                      :: caxis   = ' '
    REAL(dp), DIMENSION(:), ALLOCATABLE   :: help
    TYPE (t_geohybgrid)                   :: intgrid ! intermediate grid info
    LOGICAL                               :: lrot = .FALSE.
    INTEGER                               :: nvar, nx, ny, nn, nz
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: press4d => NULL()
    TYPE(t_chaobj_cpl)                    :: pressiobj
    TYPE(t_chaobj_cpl)                    :: pressmobj
    TYPE(t_scrip_data), POINTER           :: PSD
    INTEGER                               :: reprid = 0
    INTEGER                               :: SCRIP_ID     ! SCRIP DATA ID
    REAL(dp), DIMENSION(:,:,:), POINTER   :: vax_in  => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER   :: vax_out => NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: vdat    => NULL()
    TYPE(t_geohybgrid)                    :: xgrid
    INTEGER                               :: xsize, ysize
    INTEGER                               :: zdim
    INTEGER                               :: igridid
#endif

    IF (l_tracer_initfromrestart .and. lstart) THEN
       CALL warning_bi('Tracer init from restart file:             ', substr)
       CALL warning_bi('This is a hack! Use it very carefully!     ', substr)
       CALL warning_bi('The dimensions need to fit perfectly, i.e.,', substr)
       CALL warning_bi('the same resolution is mandatory,          ', substr)
       CALL warning_bi('and exactly the same domain!               ', substr)
       CALL warning_bi('It does only work for the GP tracer set.   ', substr)
       CALL warning_bi('It does NOT work for ICON.                 ', substr)

       CALL main_channel_read_restart('tracer_gp')

       DO ix = 1, SIZE(ti_gp)
          CALL get_channel_object_info(status, 'tracer_gp'  &
               , TRIM(ti_gp(ix)%tp%ident%fullname), lrestart_read=linit )
          CALL channel_halt(substr, status)
          ti_gp(ix)%tp%meta%linit = ti_gp(ix)%tp%meta%linit .OR. linit
       END DO
    END IF

    np = dom_curid

#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
    ldompar = .FALSE.
#ifdef PNCREGRID
    lvarpar = .TRUE.
#else
    lvarpar = .FALSE.
#endif
#endif

#if defined(COSMO) || (defined(CESM1) && defined(HOMMESE))
    ldompar = .TRUE.
    lvarpar = .FALSE.

#ifdef COSMO
    ! UPDATE PRESSURE FIELD OF BASEGRID IF REQUIRED
    pressiobj%CHA = basemodstr
    pressiobj%OBJ = 'pressi'
    pressmobj%CHA = basemodstr
    pressmobj%OBJ = 'press'
#endif

    CALL LOCATE_GEOHYBGRID(status, BASEGRID_ID(np), grid=xgrid)
    IF (status /= 0) CALL error_bi('BASEGRID NOT LOCATABLE',substr)
    ! GET PRESSURE FIELD, if base grid is height grid
    ! => ECHAM and CESM do not need these fields (xgrid%lpres = F)
    ! => COSMO and ICON require it
    IF (QDEF_NCVAR(xgrid%pressm)) THEN ! is pressm defined?
       ! GET PRESSURE FIELD
       CALL GET_CHANNEL_OBJECT(status, TRIM(pressmobj%cha),TRIM(pressmobj%OBJ) &
            , p4=press4d, dom_id=np)
       CALL CHANNEL_HALT(substr, status)
       CALL GET_CHANNEL_OBJECT_INFO(status                           &
            , TRIM(pressmobj%cha),TRIM(pressmobj%obj), reprid=reprid &
            , dom_id=np )
       CALL CHANNEL_HALT(substr, status)
       CALL GET_REPRESENTATION_INFO(status, ' ', ID=reprid, axis=caxis)
       CALL CHANNEL_HALT(substr, status)
       CALL lcase(caxis)
       CALL RGTOOL_CONVERT_DAT2PREDEFVAR(status, xgrid%pressm &
            , press4d, caxis, 'xyzn')
       CALL GEOHYBGRID_UPDATE_PRESS(status, BASEGRID_ID(np) &
               , pressm=xgrid%pressm )
       NULLIFY(press4d)
    END IF

    IF (QDEF_NCVAR(xgrid%pressi)) THEN ! is pressi defined?
       ! GET PRESSURE FIELD
       CALL GET_CHANNEL_OBJECT(status, TRIM(pressiobj%cha),TRIM(pressiobj%OBJ) &
            ,p4=press4d, dom_id=np)
       CALL CHANNEL_HALT(substr, status)
       CALL GET_CHANNEL_OBJECT_INFO(status &
            , TRIM(pressiobj%cha),TRIM(pressiobj%obj), reprid=reprid,dom_id=np)
       CALL CHANNEL_HALT(substr, status)
       CALL GET_REPRESENTATION_INFO(status, ' ', ID=reprid, axis=caxis)
       CALL CHANNEL_HALT(substr, status)
       CALL lcase(caxis)
       CALL RGTOOL_CONVERT_DAT2PREDEFVAR(status, xgrid%pressi &
            , press4d, caxis, 'xyzn')
       CALL GEOHYBGRID_UPDATE_PRESS(status, BASEGRID_ID(np) &
               , pressi=xgrid%pressi )
       NULLIFY(press4d)
    END IF
    CALL INIT_GEOHYBGRID(xgrid)
#endif
!end of #if defined(COSMO) || (defined(CESM1) && defined(HOMMESE))

    CALL start_message_bi(substr,'TRACER INITIALISATION',substr)

    ! CHECK IF PROCEDURE CAN BE SKIPPED, E.G. AFTER RESTART WITHOUT
    ! ADDING NEW TRACERS ...
    ! NOTE: lforce_init HAS BEEN TESTED ALREADY !!!
    !       HERE IS TEST OF linit SUFFICIENT !!!
    lskip = lresume ! skip only after restart
    set_loop1: DO n=1, NSETID

       CALL get_tracer_set(status, n, setname, trlist=til, ntrac=ntrac &
            , l_init=l_init)
       CALL tracer_halt(substr, status)

       ! NO EMPTY SETS
       IF (ntrac == 0) CYCLE

       ! INITIALISATION MUST BE ALLOWED
       IF (.NOT. l_init) CYCLE

       ! CHECK FOR UNINITIALIZED TRACERS
       DO
          IF (.NOT. ASSOCIATED(til)) EXIT
          ! i = til%info%ident%idx
          IF (.NOT. til%info%meta%linit) lskip = .false.
          til => til%next
       END DO

    END DO set_loop1

    IF (lskip) THEN
       CALL end_message_bi(substr,'TRACER INITIALISATION (SKIPPED)',substr)
       RETURN
    END IF

    np = dom_curid

    ! LOCATE BASEMODEL GRID
    CALL INIT_GEOHYBGRID(ogrid)
    CALL LOCATE_GEOHYBGRID(status, BASEGRID_ID(np), grid=ogrid)

    IF (status /= 0) CALL error_bi(grid_error(status), substr)

    INQUIRE(file=TRIM(modstr)//'.nml', exist=lex)  ! now tracer.nml
    IF (lex) THEN
       iunit = find_next_free_unit(100,200)
    ELSE
       IF (p_parallel_io) THEN
          CALL RGMSG(substr, RGMLW, &
               'NAMELIST FILE '''//TRIM(modstr)//'.nml'' NOT FOUND !')
          CALL RGMSG(substr, RGMLWC, &
               'NO TRACER INITIALIZATION POSSIBLE !')
       END IF
    END IF

    IF (.NOT.lex) RETURN

    ! INIT
    ALLOCATE(ni(NSETID))
    ALLOCATE(nu(NSETID))
    ALLOCATE(nr(NSETID))
    nu(:) = 0
    nr(:) = 0
    ni(:) = 0

#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
    ! ALLOCATE SPACE FOR I/O (global field)
    IF (p_parallel_io) THEN
       ALLOCATE(zin(nlon, nlev, ngl), STAT=status)
       CALL ERRMSG(substr,status,1)
    ELSE
       NULLIFY(zin)
    ENDIF
#endif

#if defined(COSMO)
    ! ALLOCATE SPACE FOR I/O (global field)
    ALLOCATE(zin(_RI_XYZ__(nproma,ngpblks,nlev)), STAT=status)
    CALL ERRMSG(substr,status,1)
#endif

    ! ALLOCATE decomposed data field (local)
    ALLOCATE(zinl(_RI_XYZ__(nproma,ngpblks,nlev)), STAT=status)
    CALL ERRMSG(substr,status,2)

    ! START REGRIDDING
    ! EXAMPLE: REGRIDDING ALL VARIABLES IN ALL NAMELISTS
    RG_NML  = NML_NEXT  ! READ NEXT NAMELIST FROM FILE
    RG_CTRL = RG_SCAN

    regrid_loop: DO ! ENDLESS DO LOOP (MUST BE TERMINATED WITH EXIT)

       RG_STATUS = RGSTAT_NULL
       CALL INIT_PIMPGRID
       CALL INIT_PARALLEL(p_pe, p_nprocs, p_all_comm)

       ! INITIALISE
       CALL INIT_GEOHYBGRID(igrid)
       tmin  = t_undef
       tc    = t_undef

       IF (ldompar .OR. lvarpar .OR. p_parallel_io) THEN

          CALL READ_CONTROL(RG_CTRL, RG_NML, RG_STATUS    &
               , var, igrid, ogrid, RGT, lint, TRIM(modstr)//'.nml' &
               , infostr                                       &
               , iounit=iunit, ldomainpar=ldompar, lvarparallel=lvarpar &
               , lpresaxis=lpres, lwork=lwork)

       END IF

       CALL p_bcast (RG_STATUS, p_io)

       IF (RG_STATUS == RGSTAT_STOP) EXIT

       outer: IF (ldompar .OR. lvarpar .OR. p_parallel_io) THEN

#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
          ! NCREGRID INTERPOLATION
          IF (lwork) &
               CALL REGRID_CONTROL(igrid, ogrid, var, ovar, RGT, lint &
               , lfirsto=(tc==tmin))
#else

          if_work: IF (lwork) THEN

             ! SCRIP INTERPOLATION
             CALL COMPLETE_GEOHYBGRID(igrid)

             CALL NEW_GEOHYBGRID(status, igridid, igrid)
             IF (status /= 0 .AND. status /= 1) &
                  CALL error_bi(grid_error(status),substr)

             ! CALCULATE SCRIP DATA AND WEIGHTS
             CALL CALC_SCRIPDATA(status, igrid, ogrid, RGT, SCRIP_ID &
                  , PSD=PSD)
             IF (status /= 0 .AND. status /= 1) &
                  CALL error_bi(grid_error(status),substr)

             IF (status == 0) THEN
                CALL CALC_SCRIP_WEIGHTS(status, PSD)
                IF (status /= 0) CALL error_bi(grid_error(status), substr)
             END IF
             NULLIFY(PSD)

             CALL SCRIP_CONTROL(status, SCRIP_ID, igrid, ogrid, RGT, lint&
                  , var, sovar,intgrid, llrgz=.TRUE., lfirsto=(tc==tmin))
             IF (status /= 0) CALL error_bi(grid_error(status), substr)

             ! make grid with vertical input grid, but horizontal
             ! SCRIP-gridded grid

             ! CONVERT VARIABLE into 4D SPACE
             CALL RGMSG(substr, RGMLI, 'CONVERT SOVAR for hori. dims')
             CALL RGTOOL_CONVERT(sovar(1), dat, intgrid,order='xyzn')
             xsize = SIZE(dat,1)
             ysize = SIZE(dat,2)
             DEALLOCATE(dat)
             NULLIFY(dat)
             ! DEFINITION OF VERTICAL AXIS: IN-GRID
             CALL RGMSG(substr, RGMLI, 'DEFINE VERTICAL InGrid Axis')
             !
             CALL CONSTRUCT_VERTICAL_AXIS(status, xsize, ysize, vax_in &
                  , igrid, lpres, SCRIP_ID, ogrid, RGT, lint)
             IF (status /= 0) &
                  CALL ERRMSG('CONSTRUCT_VERTICAL_AXIS: ' ,status,1)

             ! DEFINITION OF VERTICAL AXIS: OUT-GRID
             CALL RGMSG(substr, RGMLI, 'DEFINE VERTICAL OutGrid Axis')

             ! CHECK FOR SPECIAL CASE, vertical axis defined with hybrid
             ! coefficients without surface pressure => get surface pressure
                ! from input grid (if available)
             IF (.NOT.(QDEF_NCVAR(ogrid%pressi) &
                  .OR. QDEF_NCVAR(ogrid%pressm))) THEN
                IF (QDEF_NCVAR(ogrid%hybi)  &
                     .AND.( .NOT. QDEF_NCVAR(ogrid%ps))) THEN
                   IF (QDEF_NCVAR(igrid%ps)) THEN
                      CALL INTERPOL_GEOHYBGRID_PS(status, igrid, ogrid &
                           , SCRIP_ID)
                   ELSE
                      CALL ERRMSG(&
                      'tracer_init: WRONG INPUT FOR VERTICAL AXIS DEFINITION: '&
                           ,42,1)
                   END IF
                END IF
             END IF
             CALL CONSTRUCT_VERTICAL_AXIS(status, xsize, ysize, vax_out &
                  , ogrid, lpres)
             IF (status /= 0) &
                  CALL ERRMSG('CONSTRUCT_VERTICAL_AXIS: ' ,status,1)
             ! CHECK if both vertical axis are orientated in the same way
             ! assume axis orientation is equal for all columns, check only
             ! point (1,1)
             IF ( ( (vax_in(1,1,1)-vax_in(1,1,SIZE(vax_in,3))) * &
                  (vax_out(1,1,1)-vax_out(1,1,SIZE(vax_out,3)))) < 0) THEN
                ! V-AXIS orientation differs
                ! => AXIS ROTATION in VAX_IN and DAT required
                lrot = .TRUE.
             ELSE
                lrot = .FALSE.
             END IF
             ALLOCATE(ovar(SIZE(sovar)))
             numvar: DO nvar =1, SIZE(sovar)
                ! CONVERT VARIABLE ON INTGRID TO 3D
                ! NOTE: vertical coordinate is in intgrid 'n'
                CALL RGTOOL_CONVERT(sovar(nvar), dat, intgrid,order='xyzn')
                ! ALLOCATE MEMORY FOR VERTICAL REMAPPED FIELD
                ALLOCATE(vdat(SIZE(dat,1),SIZE(dat,2) &
                     ,SIZE(vax_out,3)-1,SIZE(dat,4)))
                vdat = 0._dp
                ! ROTATE VAX-IN and DATA IF REQUIRED
                IF (lrot) THEN
                   zdim = SIZE(VAX_IN,3)
                   ALLOCATE(help(zdim))
                   DO nx = 1, SIZE(dat,1)
                      DO ny = 1, SIZE(dat,2)
                         help(:) = VAX_IN(nx,ny,:)
                         DO nz = 1, SIZE(help)
                            VAX_IN(nx,ny,nz) = help(zdim-nz+1)
                         END DO
                         DO nn = 1, SIZE(dat,4)
                            help(1:zdim-1) = dat(nx,ny,:,nn)
                            DO nz = 1, zdim-1
                               dat(nx,ny,nz,nn) = help(zdim-1-nz+1)
                            END DO
                         END DO
                      END DO
                   END DO
                   DEALLOCATE(help)
                END IF
                DO nx = 1, SIZE(dat,1)
                   DO ny = 1, SIZE(dat,2)
                      DO nn = 1, SIZE(dat,4)
                         !
                         CALL SNREGRID(status                             &
                              , vax_in(nx,ny,:), vax_out(nx,ny,:)         &
                              , dat(nx,ny,:,nn), vdat(nx,ny,:,nn), .FALSE.)
                         IF (status /= 0) &
                              CALL ERRMSG('SNREGRID ERROR: ' ,status,25)
                      END DO
                   END DO
                END DO
                ! convert vdat to var ..
                CALL RGTOOL_CONVERT_DAT2VAR(ovar(nvar), vdat &
                     , var(nvar)%name, ogrid, 'xyzn')
                ! free memory
                DEALLOCATE(dat)
                NULLIFY(dat)
                DEALLOCATE(vdat)
                NULLIFY(vdat)
             END DO numvar
             DEALLOCATE(vax_in)
             NULLIFY(vax_in)
             DEALLOCATE(vax_out)
             NULLIFY(vax_out)

             CALL INIT_GEOHYBGRID(intgrid)
          ENDIF if_work
#endif
! end of else: #if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))

          IF (ASSOCIATED(RGT)) THEN
             DEALLOCATE(RGT) ; NULLIFY(RGT)
          ENDIF

       ENDIF outer

       do_work: DO work_loop = 0, nproc_work-1
          IF (lvarpar) THEN
             root_pe = pe_list(work_loop)
          ELSE
             root_pe = p_io
          ENDIF

          IF (ldompar) THEN
             isizev = SIZE(ovar)
          ELSE
             IF ( p_pe == root_pe ) isizev = SIZE(ovar)
             CALL p_bcast (isizev, root_pe)
          END IF

          IF ( p_pe == root_pe ) THEN
             CALL RGMSG(substr, RGMLVL, '')
             CALL RGMSG(substr, RGMLVL, 'SUBROUTINE '//TRIM(substr)//' REPORT:')
             CALL RGMSG(substr, RGMLVL, &
                  'FOUND ',isizev,' FIELD(S); CHECKING FOR TRACER ...')
          END IF

          l_proc = ldompar .OR.                             &
               ( (lvarpar .AND. p_pe == root_pe)        &
               .OR. (.NOT.lvarpar .AND. p_parallel_io) )

          species_loop: DO i = 1, isizev ! LOOP OVER SPECIES

             IF (l_proc) THEN

                CALL RGTOOL_CONVERT(ovar(i), dat, ogrid, order=repr_def_axes(_RI_XYZN_('x','y','z','n')))
                cevar = TRIM(var(i)%name)

             END IF

             IF (.NOT. ldompar) THEN
                CALL p_bcast(cevar, root_pe)
             END IF

             IF (.NOT. lvarpar) THEN
                IF (.NOT. ldompar) THEN
                   IF (p_pe == root_pe) zin(:,:,:) = dat(_RI_XYZN_(:,:,:,1))
                ENDIF
             ELSE
                ! lvarpar = T
                IF (p_pe == root_pe) THEN

                   IF(p_nprocs == 1 .OR. p_parallel_io)   THEN
                      zin(:,:,:) = dat(_RI_XYZN_(:,:,:,1))
                   ELSE
#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
                      CALL p_send(dat(_RI_XYZN_(:,:,:,1)), p_io, jt)
#endif
                   END IF
                ELSE IF (p_parallel_io .AND. p_pe /= root_pe)   THEN
                   IF(p_nprocs > 1)   THEN
#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
                      CALL p_recv(zin, root_pe, jt)
#endif
                   END IF
                END IF
             END IF

             zinl(:,:,:) = 0.0_dp
             IF (ldompar) THEN
                zinl(:,:,:) = dat(_RI_XYZN_(:,:,:,1))
             ELSE
                CALL scatter_gp(zin, zinl, dcg)
             END IF

             ! HOOK FOR INITIALISATION OF PROGNOSTIC VARIABLES
             CALL init_pvar(TRIM(cevar), zinl)

             ! TRACER HANDLING STARTS HERE ...
             set_loop2: DO n=1, NSETID

                CALL get_tracer_set(status, n, setname, ti=ti, ntrac=ntrac &
                     , l_init=l_init)
                CALL tracer_halt(substr, status)

                ! NO EMPTY SETS
                IF (ntrac == 0) CYCLE

                ! INITIALISATION MUST BE ALLOWED
                IF (.NOT. l_init) CYCLE

                CALL full2base_sub(status, TRIM(cevar), basename, subname)
                CALL tracer_halt(substr, status)
                !
                CALL get_tracer(status, setname               &
                     , TRIM(basename), TRIM(subname), idx=jt)
                !
                tracer_exists: IF (status == TR_NEXIST) THEN
                   IF (p_pe == root_pe) THEN
                      CALL RGMSG(substr, RGMLVL, &
                           '  TRACER '''//TRIM(cevar)//&
                           &''' (SET '//TRIM(setname)//') NOT DEFINED')
                   END IF
                   nu(n) = nu(n) + 1
                ELSE IF (status /= 0) THEN
                   CALL tracer_halt(substr, status)
                ELSE
                   ! skip if already initialised (from restart-file)
                   IF (ti(jt)%tp%meta%linit) THEN
                      IF ( p_pe == root_pe) THEN
                         CALL RGMSG(substr, RGMLVL, &
                              '  TRACER '''//TRIM(var(i)%name)&
                              &//''' (SET '//TRIM(setname)//&
                              &') ALREADY INITIALIZED (e.g. FROM RESTART)')
                      END IF
                      nr(n) = nr(n) + 1
                      CYCLE ! next species
                   END IF
                   !
                   IF ( p_pe == root_pe) THEN
                      CALL RGMSG(substr, RGMLVL, &
                           '  INITIALIZING TRACER '''//TRIM(var(i)%name)//&
                           &''' (SET '//TRIM(setname)//') ')
                   ENDIF

                   SELECT CASE (TRIM(setname))
                      !
                   CASE(GPTRSTR)
                      !
                      xt(_RI_XYZN_(:,:,:,jt)) = zinl(:,:,:)
                      IF (lresume) xtm1(_RI_XYZN_(:,:,:,jt)) = &
                           (1._DP - eps) * xt(_RI_XYZN_(:,:,:,jt))
                      ! set flag indicating that tracer is already initialized
                      ti(jt)%tp%meta%linit = .TRUE.
                      !
#ifdef ECHAM5
                   CASE(LGTRSTR)
                      !
!!#D attila +
                      ! distribute tracers over the processors (GP->LG)
                      zxt_lg => xt_a(_RI_XYZN_(:,1,1,jt))
                      CALL gp2lg_e5(zinl, zxt_lg)
                      IF (lresume) xtm1_a(_RI_XYZN_(:,:,:,jt)) = &
                           (1._DP - eps) * xt_a(_RI_XYZN_(:,:,:,jt))
                      ! set flag indicating that tracer is already initialized
                      ti(jt)%tp%meta%linit = .TRUE.
!!#D attila -

!!#D clams +
                   CASE(CLTRSTR)
                      ! distribute tracers over the processors (GP->CL)
                      zxt_cl => xt_c(_RI_XYZN_(:,1,1,jt))
                      CALL gp2cl_e5(zinl, zxt_cl)
! Not leap-frog for CLaMS tracers ...
!!$                   IF (lresume) xtm1_cl(_RI_XYZN_(:,:,:,jt)) = &
!!$                        (1._DP - eps) * xt_c(_RI_XYZN_(:,:,:,jt))
                      ! set flag indicating that tracer is already initialized
                      ti(jt)%tp%meta%linit = .TRUE.
!!#D clams -
#endif
                      !
                   END SELECT

                   ni(n) = ni(n) + 1

                END IF tracer_exists

             END DO set_loop2

          END DO species_loop

       END DO do_work

       IF (ASSOCIATED(var)) THEN
          DO i=1,SIZE(var)
             CALL INIT_NCVAR(var(i))
          END DO
          DEALLOCATE(var)
       END IF
       NULLIFY(var)
       IF (ASSOCIATED(dat)) THEN
          DEALLOCATE(dat, STAT=status)
          CALL ERRMSG(substr,status,4)
       END IF
       NULLIFY(dat)
       IF (ASSOCIATED(ovar)) THEN
          DO i=1, SIZE(ovar)
             CALL INIT_NCVAR(ovar(i))
          END DO
          DEALLOCATE(ovar)
          NULLIFY(ovar)
       ENDIF

       IF (ASSOCIATED(sovar)) THEN
          DO i=1, SIZE(sovar)
             CALL INIT_NCVAR(sovar(i))
          END DO
          DEALLOCATE(sovar)
          NULLIFY(sovar)
       ENDIF

       set_loop3: DO n=1, NSETID

          CALL get_tracer_set(status, n, setname &
               , ntrac=ntrac, l_init=l_init)
          CALL tracer_halt(substr, status)
          IF (ntrac == 0) CYCLE
          IF (.NOT. l_init) CYCLE

          IF (p_parallel_io) THEN
             CALL RGMSG(substr, RGMLVL, &
                  'TRACER SET : '//TRIM(setname))
             CALL RGMSG(substr, RGMLVL, &
                  '... ',nu(n), &
                  ' UNRECOGNIZED NAMES')
             CALL RGMSG(substr, RGMLVL, &
                  '... ',nr(n), &
                  ' TRACER(S) ALREADY INITIALIZED (e.g. FROM RESTART)')
             CALL RGMSG(substr, RGMLVL, &
                  '... ',ni(n),' TRACER(S) INITITALIZED')
          END IF

       END DO set_loop3

       IF (p_parallel_io) THEN
          CALL RGMSG(substr, RGMLVL, 'END OF SUBROUTINE '&
               &//TRIM(substr)//' REPORT!')
          CALL RGMSG(substr, RGMLVL, '')
       END IF

    END DO regrid_loop ! ENDLESS DO LOOP
    ! END REGRIDDING

    ! CLEAN
    CALL INIT_GEOHYBGRID(ogrid)
    CALL INIT_GEOHYBGRID(igrid)

    IF (ASSOCIATED(zin)) THEN
       DEALLOCATE(zin, STAT=status)
       CALL ERRMSG(substr,status,5)
    END IF
    NULLIFY(zin)

    IF (ASSOCIATED(zinl)) THEN
       DEALLOCATE(zinl, STAT=status)
       CALL ERRMSG(substr,status,6)
    END IF
    NULLIFY(zinl)

    DEALLOCATE(nu)
    DEALLOCATE(nr)
    DEALLOCATE(ni)

    CALL end_message_bi(modstr,'TRACER INITIALISATION',substr)
#endif
!  ... of #if defined(ECHAM5) || defined(COSMO) || defined(CESM1)

! ICON
#if defined(ICON)

    ! INITIALIZES TRACER FIELDS FROM netCDF FILES
    ! (USING IMPORT_GRID) AS DEFINED IN NAMELIST FILE tracer.nml'
    !
    ! Author: Patrick Joeckel, MPICH, Mainz, November 2002
    !         Patrick Joeckel, MPICH, Mainz, March    2004
    !         Patrick Joeckel, DLR-IPA, Oberpfaffenhofen, December 2009
    !         Astrid  Kerkweg, UNI-Mz, Mainz, August 2013
    !         Astrid  Kerkweg, UNI-BN, Bonn, September 2017

    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_bcast  &
                                      , p_pe, p_nprocs, p_all_comm
    USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                      , error_bi, info_bi
    ! ... GPTRSTR
    USE messy_main_bmluse_bi,     ONLY: p_patch, nproma
    USE messy_main_grid_def_mem_bi, ONLY: BASEGRID_ID
    USE messy_main_tools,         ONLY: find_next_free_unit, PTR_3D_ARRAY &
                                      , lcase, strcrack                   &
                                      , str2num
    USE messy_main_timer,         ONLY: lresume

    ! ... TRACER SETS
    USE messy_main_tracer,        ONLY: NSETID, STRLEN_TRSET, get_tracer_set &
                                      , t_trinfo_list, TR_NEXIST
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, STRLEN_ULONG

    ! ... IMPORT / GRID_TRAFO
    USE messy_main_import_grid,     ONLY: RG_CTRL, RG_NML              &
                                        , NML_NEXT, NML_STAY  &
                                        , RG_SCAN, RG_STATUS           &
                                        , RGSTAT_STOP, READ_CONTROL
    USE messy_main_import_grid_par, ONLY: RGSTAT_NULL, INIT_PIMPGRID   &
                                        , INIT_PARALLEL, pe_list
    USE messy_main_grid,            ONLY: t_geohybgrid, INIT_GEOHYBGRID &
                                        , grid_error, LOCATE_GEOHYBGRID &
                                        , NEW_GEOHYBGRID              &
                                        , LOCATE_GEOHYBGRID           &
                                        , INIT_GEOHYBGRID             &
                                        , GEOHYBGRID_UPDATE_PRESS
    USE messy_main_grid_netcdf,     ONLY: RGMSG, ERRMSG, RGMLVL, RGMLW     &
                                        , RGMLWC, RGMLI, INIT_NCVAR &
                                        , GRD_MAXSTRLEN, t_ncvar       &
                                        , QDEF_NCVAR
     USE messy_main_grid_tools,      ONLY: RGTOOL_CONVERT, RGTOOL_CONVERT_DAT2VAR&
                                        , RGTOOL_CONVERT_DAT2PREDEFVAR
    USE messy_main_grid_trafo,      ONLY: COMPLETE_GEOHYBGRID
    USE messy_main_grid_trafo_scrp, ONLY: t_scrip_data, SCRIP_CONTROL        &
                                        , CALC_SCRIPDATA, CALC_SCRIP_WEIGHTS &
                                        , CONSTRUCT_VERTICAL_AXIS            &
                                        , INTERPOL_GEOHYBGRID_PS
    USE messy_main_grid_trafo_nrgd_base, ONLY: SNREGRID

    USE messy_main_channel_error_bi,     ONLY: CHANNEL_HALT
    USE messy_main_channel_bi,           ONLY: n_dom
    USE messy_main_channel,              ONLY: GET_CHANNEL_OBJECT      &
                                             , t_chaobj_cpl            &
                                             , GET_CHANNEL_OBJECT_INFO
    USE messy_main_channel_repr,         ONLY: GET_REPRESENTATION_INFO
    USE messy_main_tools,                ONLY: split_name_domain

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, NULL, SIZE, TRIM

    TYPE t_ncvar1darray
       TYPE (t_ncvar), DIMENSION(:), POINTER :: ncv => NULL()
    END type t_ncvar1darray

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'tracer_init'
    INTEGER                      :: n       ! set counter
    INTEGER                      :: ntrac   ! number of tracers
    CHARACTER(LEN=STRLEN_TRSET)  :: setname ! tracer set name
    LOGICAL                      :: l_init  ! initialize this set ?
    TYPE(t_trinfo_list), POINTER :: til => NULL()
    TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti  => NULL()
    !
    INTEGER                      :: status  ! error status
    LOGICAL                      :: lskip   ! skip, if no uninitialized tracers
    INTEGER                      :: iunit   ! fortran unit for input file
    LOGICAL                      :: lex     ! file exists ?
    !
    !
    ! decomposed data field (local)
    TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER :: zinl  => NULL()
    ! REAL(DP), DIMENSION(:,:,:), POINTER :: zinl => NULL()
    !
    ! GRID
    TYPE (t_ncvar), DIMENSION(:), POINTER :: var => NULL() ! list of variables
    TYPE (t_geohybgrid)                   :: grid  ! struct with grid info
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: dat => NULL()
    CHARACTER(len=2*STRLEN_MEDIUM + 1)    :: cevar ! object of evar
    CHARACTER(len=STRLEN_MEDIUM)          :: basename
    CHARACTER(len=STRLEN_MEDIUM)          :: subname
    INTEGER                               :: isizev ! SIZE(evar)
    INTEGER                               :: i      ! species counter
    !
    INTEGER, DIMENSION(:), ALLOCATABLE    :: nu, nr, ni ! counter
    INTEGER                               :: jt         ! tracer index in set
    INTEGER                               :: work_loop

    LOGICAL,  SAVE                        :: lint    ! use input time ?
    INTEGER, DIMENSION(:), POINTER        :: RGT => NULL() ! regridding type
    TYPE (t_ncvar1darray), DIMENSION(:), POINTER :: ovar  => NULL() ! list of variables
    TYPE (t_ncvar), DIMENSION(:), POINTER :: sovar => NULL() ! list of variables
    TYPE (t_geohybgrid)                   :: igrid   ! input grid info
    TYPE (t_geohybgrid), DIMENSION(:), POINTER :: ogrid   ! output grid info
    TYPE (t_geohybgrid)                   :: intgrid ! intermediate grid info
    INTEGER                               :: count = 0 ! procedure counter
    TYPE(t_scrip_data), POINTER           :: PSD
    INTEGER                               :: igridid
    INTEGER                               :: SCRIP_ID     ! SCRIP DATA ID
    INTEGER                               :: tc, tmin
    INTEGER,           PARAMETER          :: t_undef = -99

    LOGICAL                               :: lpres   = .FALSE.
    LOGICAL                               :: lwork   = .FALSE.
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: vdat => NULL()
    INTEGER                               :: nvar, nx, ny, nn, nz
    REAL(dp), DIMENSION(:,:,:), POINTER   :: vax_in  => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER   :: vax_out => NULL()
    INTEGER                               :: xsize, ysize
    LOGICAL                               :: lrot = .FALSE.
    INTEGER                               :: zdim
    REAL(dp), DIMENSION(:), ALLOCATABLE   :: help
    CHARACTER(LEN=4)            :: caxis   = ' '
    TYPE(t_geohybgrid)          :: xgrid
    INTEGER                     :: reprid = 0
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: press4d => NULL()
    TYPE(t_chaobj_cpl)          :: pressiobj
    TYPE(t_chaobj_cpl)          :: pressmobj
    INTEGER                     :: np
    CHARACTER(LEN=STRLEN_TRSET), DIMENSION(:), POINTER :: arrstr => NULL()
    INTEGER                      :: ix
    INTEGER                      :: num
    CHARACTER(LEN=STRLEN_TRSET)  :: setstring = ''
    CHARACTER(LEN=STRLEN_ULONG)  :: infostr = ''
    CHARACTER(LEN=STRLEN_TRSET)  :: str     = ''
    LOGICAL, DIMENSION(:), ALLOCATABLE :: linitpatchX

    CALL start_message_bi(substr,'TRACER INITIALISATION',substr)

    ! CHECK IF PROCEDURE CAN BE SKIPPED, E.G. AFTER RESTART WITHOUT
    ! ADDING NEW TRACERS ...
    ! NOTE: lforce_init HAS BEEN TESTED ALREADY !!!
    !       HERE IS TEST OF linit SUFFICIENT !!!
    lskip = lresume ! skip only after restart
    set_loop1: DO n=1, NSETID

       CALL get_tracer_set(status, n, setname, trlist=til, ntrac=ntrac &
            , l_init=l_init)
       CALL tracer_halt(substr, status)

       ! NO EMPTY SETS
       IF (ntrac == 0) CYCLE

       ! INITIALISATION MUST BE ALLOWED
       IF (.NOT. l_init) CYCLE

       ! CHECK FOR UNINITIALIZED TRACERS
       DO
          IF (.NOT. ASSOCIATED(til)) EXIT
          ! i = til%info%ident%idx
          IF (.NOT. til%info%meta%linit) lskip = .false.
          til => til%next
       END DO

    END DO set_loop1

    IF (lskip) THEN
       CALL end_message_bi(substr,'TRACER INITIALISATION (SKIPPED)',substr)
       RETURN
    END IF

    INQUIRE(file=TRIM(modstr)//'.nml', exist=lex)  ! now tracer.nml
    IF (lex) THEN
       iunit = find_next_free_unit(100,200)
    ELSE
       IF (p_parallel_io) THEN
          CALL RGMSG(substr, RGMLW, &
               'NAMELIST FILE '''//TRIM(modstr)//'.nml'' NOT FOUND !')
          CALL RGMSG(substr, RGMLWC, &
               'NO TRACER INITIALIZATION POSSIBLE !')
       END IF
    END IF

    IF (.NOT.lex) RETURN

    ! INIT
    ALLOCATE(ni(NSETID))
    ALLOCATE(nu(NSETID))
    ALLOCATE(nr(NSETID))
    nu(:) = 0
    nr(:) = 0
    ni(:) = 0

    ! UPDATE PRESSURE FIELD OF BASEGRID IF REQUIRED
    pressmobj%CHA = 'nh_state_diag'
    pressmobj%OBJ = 'pres'
    pressiobj%CHA = 'nh_state_diag'
    pressiobj%OBJ = 'pres_ifc'

    DO np = 1, n_dom
       CALL LOCATE_GEOHYBGRID(status, BASEGRID_ID(np), grid=xgrid)
       IF (status /= 0) CALL error_bi('BASEGRID NOT LOCATABLE',substr)
       ! GET PRESSURE FIELD, if base grid is height grid
       ! => ECHAM and CESM do not need these fields (xgrid%lpres = F)
       ! => COSMO and ICON require it
       IF (QDEF_NCVAR(xgrid%pressm)) THEN ! is pressm defined?
          ! GET PRESSURE FIELD
          CALL GET_CHANNEL_OBJECT(status &
               , TRIM(pressmobj%cha),TRIM(pressmobj%OBJ) &
               ,p4=press4d, dom_id=np)
          CALL CHANNEL_HALT(substr, status)
          CALL GET_CHANNEL_OBJECT_INFO(status                          &
               ,TRIM(pressmobj%cha),TRIM(pressmobj%obj), reprid=reprid &
               ,dom_id=np )
          CALL CHANNEL_HALT(substr, status)
          CALL GET_REPRESENTATION_INFO(status, ' ', ID=reprid, axis=caxis)
          CALL CHANNEL_HALT(substr, status)
          CALL lcase(caxis)
          CALL RGTOOL_CONVERT_DAT2PREDEFVAR(status, xgrid%pressm &
               , press4d, caxis, 'xyzn')
          CALL GEOHYBGRID_UPDATE_PRESS(status, BASEGRID_ID(np) &
               , pressm=xgrid%pressm )
          NULLIFY(press4d)
       END IF
       IF (QDEF_NCVAR(xgrid%pressi)) THEN ! is pressi defined?
          ! GET PRESSURE FIELD
          CALL GET_CHANNEL_OBJECT(status                 &
               , TRIM(pressiobj%cha),TRIM(pressiobj%OBJ) &
               , p4=press4d, dom_id=np)
          CALL CHANNEL_HALT(substr, status)
          CALL GET_CHANNEL_OBJECT_INFO(status            &
               ,TRIM(pressiobj%cha),TRIM(pressiobj%obj)  &
               , reprid=reprid,dom_id=np )
          CALL CHANNEL_HALT(substr, status)
          CALL GET_REPRESENTATION_INFO(status, ' ', ID=reprid, axis=caxis)
          CALL CHANNEL_HALT(substr, status)
          CALL lcase(caxis)
          CALL RGTOOL_CONVERT_DAT2PREDEFVAR(status, xgrid%pressi &
               , press4d, caxis, 'xyzn')
          CALL GEOHYBGRID_UPDATE_PRESS(status, BASEGRID_ID(np) &
               , pressi=xgrid%pressi )
          NULLIFY(press4d)
       END IF
       CALL INIT_GEOHYBGRID(xgrid)
    END DO

    ! LOCATE BASEMODEL GRIDs
    ALLOCATE(ogrid(n_dom))
    DO np = 1, n_dom
       CALL INIT_GEOHYBGRID(ogrid(np))
       CALL LOCATE_GEOHYBGRID(status, BASEGRID_ID(np), grid=ogrid(np))
       IF (status /= 0) CALL error_bi(grid_error(status), substr)
    END DO

    ALLOCATE(zinl(n_dom))
    DO np= 1,n_dom
       ! ALLOCATE decomposed data field (local)
       ALLOCATE(zinl(np)%ptr(nproma,p_patch(np)%nlev,p_patch(np)%nblks_c)&
            , STAT=status)
       CALL ERRMSG(substr,status,2)
    END DO

    ALLOCATE(ovar(n_dom))
    DO np = 1,n_dom
       NULLIFY(ovar(np)%ncv)
    END DO
    ALLOCATE(linitpatchX(n_dom))

    ! START REGRIDDING
    ! EXAMPLE: REGRIDDING ALL VARIABLES IN ALL NAMELISTS
    RG_NML  = NML_NEXT  ! READ NEXT NAMELIST FROM FILE
    RG_CTRL = RG_SCAN

    regrid_loop: DO ! ENDLESS DO LOOP (MUST BE TERMINATED WITH EXIT)

       RG_STATUS = RGSTAT_NULL
       CALL INIT_PIMPGRID
       CALL INIT_PARALLEL(p_pe, p_nprocs, p_all_comm)

       ! INITIALISE
       CALL INIT_GEOHYBGRID(igrid)
       tmin = t_undef
       tc   = t_undef

       CALL READ_CONTROL(RG_CTRL, RG_NML, RG_STATUS                  &
            , var, igrid, ogrid(1), RGT, lint, TRIM(modstr)//'.nml' &
            , infostr                                                &
            , iounit=iunit, ldomainpar=.TRUE., lvarparallel=.FALSE.  &
            , lpresaxis=lpres, lwork=lwork)

       CALL p_bcast (RG_STATUS, p_io)

       IF (RG_STATUS == RGSTAT_STOP) EXIT

       IF (TRIM(ADJUSTL(infostr)) == '') THEN
       ! EMPTY infostr => init all patches
          linitpatchX(:) = .TRUE.
       ELSE
          linitpatchX(:) = .FALSE.
          CALL strcrack(infostr, ',',arrstr,num)
          DO ix = 1, num
             str = TRIM(ADJUSTL(arrstr(ix)))
             IF (str(1:2) == '#D') THEN
                CALL str2num(str(3:), np, status)
                IF (np > n_dom) THEN
                   CALL info_bi('SKIPPING INIT REQUEST FOR DOMAIN '//TRIM(str)&
                        ,substr)
                   CYCLE
                END IF
                IF (status == 0) THEN
                   linitpatchX(np) = .TRUE.
                ELSE
                   CALL error_bi(&
                     'Error in info string from &regrid nml: '//str &
                        , substr)
                END IF
             ELSE
                CALL error_bi('UNKOWN info string in &regrid nml: '//str &
                     , substr)
             END IF
          END DO
          DEALLOCATE(arrstr)
       END IF

       dom_loop: DO np = 1, n_dom

          IF (.NOT. linitpatchX(np)) CYCLE

          RG_NML = NML_STAY
          CALL INIT_GEOHYBGRID(igrid)

          CALL READ_CONTROL(RG_CTRL, RG_NML, RG_STATUS                  &
               , var, igrid, ogrid(np), RGT, lint, TRIM(modstr)//'.nml' &
               , infostr                                                &
               , iounit=iunit, ldomainpar=.TRUE., lvarparallel=.FALSE.  &
               , lpresaxis=lpres, lwork=lwork)

          CALL p_bcast (RG_STATUS, p_io)

          IF (RG_STATUS == RGSTAT_STOP) EXIT

          ! SCRIP INTERPOLATION
          CALL COMPLETE_GEOHYBGRID(igrid)

          CALL NEW_GEOHYBGRID(status, igridid, igrid)
          IF (status /= 0 .AND. status /= 1) &
               CALL error_bi(grid_error(status),substr)

          ! CALCULATE SCRIP DATA AND WEIGHTS
          CALL CALC_SCRIPDATA(status, igrid, ogrid(np), RGT, SCRIP_ID, PSD=PSD)
          IF (status /= 0 .AND. status /= 1) &
               CALL error_bi(grid_error(status),substr)

          IF (status == 0) THEN
             CALL CALC_SCRIP_WEIGHTS(status, PSD)
             IF (status /= 0) CALL error_bi(grid_error(status), substr)
          END IF
          NULLIFY(PSD)

          CALL SCRIP_CONTROL(status, SCRIP_ID, igrid, ogrid(np), RGT, lint&
               , var, sovar,intgrid, llrgz=.TRUE., lfirsto=(tc==tmin))
          IF (status /= 0) CALL error_bi(grid_error(status), substr)

          !  grid with vertical input grid, but horizontal SCRIP-gridded grid

          ! CONVERT VARIABLE into 4D SPACE
          CALL RGMSG(substr, RGMLI, 'CONVERT SOVAR for hori. dims')
          CALL RGTOOL_CONVERT(sovar(1), dat, intgrid,order='xyzn')
          xsize = SIZE(dat,1)
          ysize = SIZE(dat,2)
          DEALLOCATE(dat)
          NULLIFY(dat)

          ! a DEFINITION OF VERTICAL AXIS: IN-GRID
          CALL RGMSG(substr, RGMLI, 'DEFINE VERTICAL InGrid Axis')
          !
          CALL CONSTRUCT_VERTICAL_AXIS(status, xsize, ysize, vax_in &
               , igrid, lpres, SCRIP_ID, ogrid(np), RGT, lint)
          IF (status /= 0) &
               CALL ERRMSG('CONSTRUCT_VERTICAL_AXIS: ' ,status,1)

          ! a DEFINITION OF VERTICAL AXIS: OUT-GRID
          CALL RGMSG(substr, RGMLI, 'DEFINE VERTICAL OutGrid Axis')
          !
          ! CHECK FOR SPECIAL CASE, vertical axis defined with hybrid
          ! coefficients without surface pressure => get surface pressure
          ! from input grid (if available)
          IF (.NOT.(QDEF_NCVAR(ogrid(np)%pressi) &
               .OR. QDEF_NCVAR(ogrid(np)%pressm))) THEN
             IF (QDEF_NCVAR(ogrid(np)%hybi)  &
                  .AND.( .NOT. QDEF_NCVAR(ogrid(np)%ps))) THEN
                IF (QDEF_NCVAR(igrid%ps)) THEN
                   CALL INTERPOL_GEOHYBGRID_PS(status, igrid, ogrid(np) &
                           , SCRIP_ID)
                ELSE
                   CALL ERRMSG(&
                    'tracer_init: WRONG INPUT FOR VERTICAL AXIS DEFINITION: '&
                           ,42,1)
                END IF
             END IF
          END IF

          CALL CONSTRUCT_VERTICAL_AXIS(status, xsize, ysize, vax_out &
               , ogrid(np), lpres)
          IF (status /= 0) CALL ERRMSG('CONSTRUCT_VERTICAL_AXIS: ' ,status,1)

          ! CHECK if both vertical axes are oriented in the same way;
          ! assume axis orientation is equal for all columns, check only
          ! point (1,1)
          IF ( ( (vax_in(1,1,1)-vax_in(1,1,SIZE(vax_in,3))) * &
               (vax_out(1,1,1)-vax_out(1,1,SIZE(vax_out,3)))) < 0) THEN
             ! V-AXIS orientation differs
             ! => AXIS ROTATION in VAX_IN and DAT required
             lrot = .TRUE.
          ELSE
             lrot = .FALSE.
          END IF

          ALLOCATE(ovar(np)%ncv(SIZE(sovar)))
          numvar: DO nvar =1, SIZE(sovar)
             ! CONVERT VARIABLE ON INTGRID TO 3D
             ! NOTE: vertical coordinate is in intgrid 'n'
             CALL RGTOOL_CONVERT(sovar(nvar), dat, intgrid,order='xyzn')
             ! ALLOCATE MEMORY FOR VERTICAL REMAPPED FIELD
             ALLOCATE(vdat(SIZE(dat,1),SIZE(dat,2) &
                  ,SIZE(vax_out,3)-1,SIZE(dat,3)))
             vdat = 0._dp
             ! ROTATE VAX-IN and DATA IF REQUIRED
             IF (lrot) THEN
                zdim = SIZE(VAX_IN,3)
                ALLOCATE(help(zdim))
                DO nx = 1, SIZE(dat,1)
                   DO ny = 1, SIZE(dat,2)
                      help(:) = VAX_IN(nx,ny,:)
                      DO nz = 1, SIZE(help)
                         VAX_IN(nx,ny,nz) = help(zdim-nz+1)
                      END DO
                      DO nn = 1, SIZE(dat,3)
                         help(1:zdim-1) = dat(nx,ny,nn,:)
                         DO nz = 1, zdim-1
                            dat(nx,ny,nn,nz) = help(zdim-1-nz+1)
                         END DO
                      END DO
                   END DO
                END DO
                DEALLOCATE(help)
             END IF
             DO nx = 1, SIZE(dat,1)
                DO ny = 1, SIZE(dat,2)
                   DO nn = 1, SIZE(dat,3)
                      !
                      CALL SNREGRID(status                             &
                           , vax_in(nx,ny,:), vax_out(nx,ny,:)         &
                           , dat(nx,ny,nn,:), vdat(nx,ny,:,nn), .FALSE.)
                      IF (status /= 0) &
                           CALL ERRMSG('SNREGRID ERROR: ' ,status,25)
                   END DO
                END DO
             END DO
             ! convert vdat to var ..
             CALL RGTOOL_CONVERT_DAT2VAR(ovar(np)%ncv(nvar), vdat &
                  , var(nvar)%name, ogrid(np), 'xyzn')
             ! free memory
             DEALLOCATE(dat)
             NULLIFY(dat)
             DEALLOCATE(vdat)
             NULLIFY(vdat)
          END DO numvar
          DEALLOCATE(vax_in)
          NULLIFY(vax_in)
          DEALLOCATE(vax_out)
          NULLIFY(vax_out)

          CALL INIT_GEOHYBGRID(intgrid)
          CALL INIT_GEOHYBGRID(igrid)

       END DO dom_loop
       RG_NML = NML_NEXT

       IF ( p_parallel_io ) THEN
          CALL RGMSG(substr, RGMLVL, '')
          CALL RGMSG(substr, RGMLVL, 'SUBROUTINE '//TRIM(substr)//' REPORT:')
       END IF
       ! TRACER HANDLING STARTS HERE ...
       set_loop2: DO n=1, NSETID

          CALL get_tracer_set(status, n, setname, ti=ti, ntrac=ntrac &
!!$            , xt=pxt, xtte=pxtte, xtm1=pxtm1 &  ! qqq
               , l_init=l_init )
          CALL tracer_halt(substr, status)

          ! NO EMPTY SETS
          IF (ntrac == 0) CYCLE

          ! INITIALISATION MUST BE ALLOWED
          IF (.NOT. l_init) CYCLE

          ! GET PATCH ID from setname
          CALL split_name_domain(status, setname, setstring, np)

          ! SKIP IF INFO in &regrid prevents this PATCH from initialisation
          IF (.NOT. ASSOCIATED(ovar(np)%ncv)) CYCLE

          isizev = SIZE(ovar(np)%ncv) ! rescue here

          IF ( p_parallel_io ) THEN
             CALL RGMSG(substr, RGMLVL, &
                  'FOUND ',isizev,' FIELD(S); CHECKING FOR TRACER ...')
             CALL RGMSG(substr, RGMLVL, &
                  '... in tracer set '//setname)
          END IF

          ! SET POINTERS TO TRACER MEMORY TO CURRENT PATCH
          CALL main_tracer_set_domain(np)

          species_loop: DO i = 1, isizev ! LOOP OVER SPECIES

             CALL RGTOOL_CONVERT(ovar(np)%ncv(i), dat, ogrid(np), order='xzny')
             cevar = TRIM(var(i)%name)
             zinl(np)%ptr(:,:,:) = dat(:,:,1,:)

!qqq
             ! HOOK FOR INITIALISATION OF PROGNOSTIC VARIABLES
             ! XXX prog. Variable initialisation not yet implemented in ICON
             !CALL init_pvar(TRIM(cevar), zinl(np)%ptr)

             CALL full2base_sub(status, TRIM(cevar), basename, subname)
             CALL tracer_halt(substr, status)
             !
             CALL get_tracer(status, setname               &
                  , TRIM(basename), TRIM(subname), idx=jt)
             !
             tracer_exists: IF (status == TR_NEXIST) THEN
                IF (p_parallel_io) THEN
                   CALL RGMSG(substr, RGMLVL, &
                        '  TRACER '''//TRIM(cevar)//&
                        &''' (SET '//TRIM(setname)//') NOT DEFINED')
                END IF
                nu(n) = nu(n) + 1
             ELSE IF (status /= 0) THEN
                CALL tracer_halt(substr, status)
             ELSE
                ! skip if already initialised (from restart-file)
                IF (ti(jt)%tp%meta%linit) THEN
                   IF ( p_parallel_io) THEN
                      CALL RGMSG(substr, RGMLVL, &
                           '  TRACER '''//TRIM(var(i)%name)&
                           &//''' (SET '//TRIM(setname)//&
                           &') ALREADY INITIALIZED (e.g. FROM RESTART)')
                   END IF
                   nr(n) = nr(n) + 1
                   CYCLE ! next species
                END IF
                !
                IF ( p_parallel_io) THEN
                   CALL RGMSG(substr, RGMLVL, &
                        '  INITIALIZING TRACER '''//TRIM(var(i)%name)//&
                        &''' (SET '//TRIM(setname)//') ')
                ENDIF

                SELECT CASE (setstring(1:2))  !'gp'
                   !
                CASE('gp')
                   !
                   xt(_RI_XYZN_(:,:,:,jt)) = zinl(np)%ptr(:,:,:)
!qqq test                   IF (lresume) xtm1(_RI_XYZN_(:,:,:,jt)) = &
                   xtm1(_RI_XYZN_(:,:,:,jt)) = &
                        xt(_RI_XYZN_(:,:,:,jt))
                   ! set flag indicating that tracer is already initialized
                   ti(jt)%tp%meta%linit = .TRUE.
                   !
                   !
                END SELECT

                ni(n) = ni(n) + 1

             END IF tracer_exists

          END DO species_loop
       END DO set_loop2

       IF (ASSOCIATED(var)) THEN
          DO i=1,SIZE(var)
             CALL INIT_NCVAR(var(i))
          END DO
          DEALLOCATE(var)
       END IF
       NULLIFY(var)
       IF (ASSOCIATED(dat)) THEN
          DEALLOCATE(dat, STAT=status)
          CALL ERRMSG(substr,status,4)
       END IF
       NULLIFY(dat)

       IF (ASSOCIATED(ovar)) THEN
          DO np = 1, n_dom
             IF (ASSOCIATED(ovar(np)%ncv)) THEN
                DO i=1, SIZE(ovar(np)%ncv)
                   CALL INIT_NCVAR(ovar(np)%ncv(i))
                END DO
                DEALLOCATE(ovar(np)%ncv)
                NULLIFY(ovar(np)%ncv)
             END IF
          END DO
       END IF
       IF (ASSOCIATED(sovar)) THEN
          DO i=1, SIZE(sovar)
             CALL INIT_NCVAR(sovar(i))
          END DO
          DEALLOCATE(sovar)
          NULLIFY(sovar)
       ENDIF

       set_loop3: DO n=1, NSETID

          CALL get_tracer_set(status, n, setname &
               , ntrac=ntrac, l_init=l_init)
          CALL tracer_halt(substr, status)

          IF (ntrac == 0) CYCLE
          IF (.NOT. l_init) CYCLE

          IF (p_parallel_io) THEN
             CALL RGMSG(substr, RGMLVL, &
                  'TRACER SET : '//TRIM(setname))
             CALL RGMSG(substr, RGMLVL, &
                  '... ',nu(n), &
                  ' UNRECOGNIZED NAMES')
             CALL RGMSG(substr, RGMLVL, &
                  '... ',nr(n), &
                  ' TRACER(S) ALREADY INITIALIZED (e.g. FROM RESTART)')
             CALL RGMSG(substr, RGMLVL, &
                  '... ',ni(n),' TRACER(S) INITITALIZED')
          END IF

       END DO set_loop3

       IF (p_parallel_io) THEN
          CALL RGMSG(substr, RGMLVL, 'END OF SUBROUTINE '&
               &//TRIM(substr)//' REPORT!')
          CALL RGMSG(substr, RGMLVL, '')
       END IF

       IF (ASSOCIATED(RGT)) THEN
          DEALLOCATE(RGT) ; NULLIFY(RGT)
       ENDIF

    END DO regrid_loop ! ENDLESS DO LOOP
    ! END REGRIDDING

    ! CLEAN
    DEALLOCATE(ovar)
    NULLIFY(ovar)

    DO np = 1, n_dom
       CALL INIT_GEOHYBGRID(ogrid(np))
    END DO
    DEALLOCATE(ogrid)
    NULLIFY(ogrid)

    DO np = 1, n_dom
       IF (ASSOCIATED(zinl(np)%ptr)) THEN
          DEALLOCATE(zinl(np)%ptr, STAT=status)
          CALL ERRMSG(substr,status,6)
       END IF
       NULLIFY(zinl(np)%ptr)
    END DO
    DEALLOCATE(zinl)
    NULLIFY(zinl)

    DEALLOCATE(nu)
    DEALLOCATE(nr)
    DEALLOCATE(ni)

    DEALLOCATE(linitpatchX)

    CALL end_message_bi(modstr,'TRACER INITIALISATION',substr)

#endif
! of #ifdef ICON

  END SUBROUTINE TRACER_INIT

  ! ------------------------------------------------------------------

#if ! defined(ICON)
  ! -------------------------------------------------------------------
  SUBROUTINE setup_tracer_set_gp

    ! SETUP MEMORY FOR TRACERS IN GRIDPOINT REPRESENTATION

    USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev, ngpblks
    USE messy_main_data_bi,         ONLY: l2tls
#ifdef COSMOv5s5
    USE messy_main_grid_def_mem_bi, ONLY: nproma_cosmo
#endif

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'setup_tracer_set_gp'
    INTEGER                        :: status
    INTEGER, DIMENSION(3)          :: dims
#if defined(COSMOv5s5)
    INTEGER                        :: dimblck
    INTEGER                        :: ntblck
#endif

#if defined(ECHAM5)
    dims(1) = nproma
    dims(2) = nlev
    dims(3) = ngpblks

    ! 4 time levels: xt, xtte, xtm1, xtf (the latter is 'extended')
    ! first 3 time levels are 'standard' (t, tendency, t-1) -> .TRUE.
    ! tracers can be initialised -> .TRUE.
    CALL setup_tracer_set(status, GPTRSTR, dims, 4, .TRUE., .TRUE.)
    CALL tracer_halt(substr, status)
#endif
#if defined(COSMO) || defined(BLANK) || defined(MBM_CLAMS) || defined(VERTICO) || defined(MBM_MPIOM) || defined(MESSYDWARF)
    dims(:) = (/ _RI_XYZ__(nproma,ngpblks,nlev) /)
    ! 5/6 levels: xt, xtte, xtm1, xtf (the latter is 'extended')
    ! first 3 time levels are 'standard' (t, tendency, t-1) -> .TRUE.
#ifdef COSMOv5s5
    dimblck = nproma_cosmo
    ntblck  = 2
#endif

    ! A) for leapfrog (.not. l2tls)
    ! pxtf contain the memory for xtf and the 2 timelevels of the
    ! boundary data (xt_bd), with xt_bd level 1,2 and pxtf level 3
    ! B) RUNGE_KUTTA
    ! pxtf contain the memory only for 2 timelevels of the boundary data (xt_bd)

    IF (.NOT. l2tls) THEN
       ! tracers can be initialised -> .TRUE.
       CALL setup_tracer_set(status, GPTRSTR, dims, 6, .TRUE., .TRUE. &
#ifdef COSMOv5s5
            , dimblck=dimblck, ntblck=ntblck &
#endif
            )
       CALL tracer_halt(substr, status)
    ELSE
       ! tracers can be initialised -> .TRUE.
       CALL setup_tracer_set(status, GPTRSTR, dims, 5, .TRUE., .TRUE. &
#ifdef COSMOv5s5
            , dimblck=dimblck, ntblck=ntblck &
#endif
            )
       CALL tracer_halt(substr, status)
    ENDIF
#endif
#if defined(CESM1)
    dims(1) = nproma
    dims(2) = nlev
    dims(3) = ngpblks

    ! 3 time levels: xt, xtte, xtm1
    ! first 3 time levels are 'standard' (t, tendency, t-1) -> .TRUE.
    ! tracers can be initialiszed -> .TRUE.
    CALL setup_tracer_set(status, GPTRSTR, dims, 3, .TRUE., .TRUE.)
    CALL tracer_halt(substr, status)
#endif
    ! SETUP POINTERS TO GRIDPOINT TRACER MEMORY AND TRACER INFORMATION
    CALL get_tracer_set(status, GPTRSTR, trlist_gp, ti_gp, ntrac_gp &
         , xt=pxt, xtte=pxtte, xtm1=pxtm1, xmem=pxtf                &
#ifdef COSMO
         , xt_tl=xt_tl &
#ifdef COSMOv5s5
         , ntracblck=ntracblck_gp &
         , xtblck=xtblck, xtteblck=xtteblck    &
#endif
#endif
         )
    CALL tracer_halt(substr, status)

  END SUBROUTINE setup_tracer_set_gp
  ! -------------------------------------------------------------------

#else
! ICON

  ! -------------------------------------------------------------------
  SUBROUTINE setup_tracer_set_gp

    ! SETUP MEMORY FOR TRACERS IN GRIDPOINT REPRESENTATION

    USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev, ngpblks
    USE messy_main_channel_mem,     ONLY: dom_current

    IMPLICIT NONE

    INTRINSIC :: NULL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'setup_tracer_set_gp'
    INTEGER                        :: status
    INTEGER, DIMENSION(3)          :: dims

    dims(1) = nproma
    dims(2) = nlev
    dims(3) = ngpblks

!qqq role of time levels needs to be checked, (how) are they "rotated"?

    ! 4 time levels: xt, xtte, xtm1, xtf (the latter is 'extended')
    ! first 3 time levels are 'standard' (t, tendency, t-1) -> .TRUE.
    ! tracers can be initialised -> .TRUE.
    CALL setup_tracer_set(status, L_GPTRSTR(dom_current), dims &
!!$      , nsav1(dom_current), .TRUE., .TRUE.)
         , 3, .TRUE., .TRUE.)
    CALL tracer_halt(substr, status)

    ! SETUP POINTERS TO GRIDPOINT TRACER MEMORY AND TRACER INFORMATION
    CALL get_tracer_set(status, L_GPTRSTR(dom_current)  &
         , trlist_gp, ti_gp, ntrac_gp                     &
         , xt=pxt, xtte=pxtte, xtm1=pxtm1)
    CALL tracer_halt(substr, status)

  END SUBROUTINE setup_tracer_set_gp
  ! -------------------------------------------------------------------
#endif

#ifdef ECHAM5
  ! -------------------------------------------------------------------
  SUBROUTINE setup_tracer_set_cl

      USE messy_main_tracer_mem_bi,   ONLY: dnparts_max
      USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev, ngpblks

    ! SETUP MEMORY FOR TRACERS IN CLaMS REPRESENTATION

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'setup_tracer_set_cl'
    INTEGER                        :: status
    INTEGER, DIMENSION(3)          :: dims_cl
    INTEGER, DIMENSION(3)          :: dims_clgp

    dims_cl(1) = dnparts_max
    dims_cl(2) = 1
    dims_cl(3) = 1

    ! 3 time levels: xt_c, xtte_c, xtm1_c
    ! These 3 time levels are 'standard' (t, tendency, t-1) -> .TRUE.
    ! tracers can be initialiszed -> .TRUE.
    CALL setup_tracer_set(status, CLTRSTR, dims_cl, 2, .FALSE., .TRUE.)
    CALL tracer_halt(substr, status)

    ! SETUP POINTERS TO CLaMS TRACER MEMORY AND TRACER INFORMATION
    CALL get_tracer_set(status, CLTRSTR, trlist_cl, ti_cl, ntrac_cl &
         , xt=pxt_c, xmem=pxmem_c) !, xtte=pxtte_c, xtm1=pxtm1_c)
    CALL tracer_halt(substr, status)

    pxtm1_c => pxt_c(:,:,:,:,:)      ! special for CLaMS
    pxtte_c => pxmem_c(:,:,:,:,1:1)  ! ext. memory for tendency

    conversion_cl2gp: IF (l_conv_cl2gp) THEN

       CALL copy_tracer_set(status, CLTRSTR, CLGPTRSTR )
       CALL tracer_halt(substr, status)

       dims_clgp(1) = nproma
       dims_clgp(2) = nlev
       dims_clgp(3) = ngpblks

       ! only one time level, no standard, tracers shall not be externally
       ! initialized
       CALL setup_tracer_set(status, CLGPTRSTR, dims_clgp, 1, .FALSE., .FALSE.)
       CALL tracer_halt(substr, status)

       CALL get_tracer_set(status, CLGPTRSTR, trlist_clgp, ti_clgp, ntrac_clgp &
            , xt=pxt_clgp)!, xtte=pxtte_clgp, xtm1=pxtm1_clgp)
       CALL tracer_halt(substr, status)

    END IF conversion_cl2gp

  END SUBROUTINE setup_tracer_set_cl
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE setup_tracer_set_lg

    ! SETUP MEMORY FOR TRACERS IN LAGRANGIAN REPRESENTATION

    USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev, ngpblks

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, NULL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'setup_tracer_set_lg'
    INTEGER                        :: status
    INTEGER, DIMENSION(3)          :: dims
    TYPE(t_trinfo_list), POINTER   :: til => NULL()

    dims(1) = NCELL
    dims(2) = 1
    dims(3) = 1

    ! 4 time levels: xt_a, xtte_a, xtm1_a, xtf_a (the latter is 'extended')
    ! first 3 time levels are 'standard' (t, tendency, t-1) -> .TRUE.
    ! tracers can be initialiszed -> .TRUE.
    CALL setup_tracer_set(status, LGTRSTR, dims, 4, .TRUE., .TRUE.)
    CALL tracer_halt(substr, status)

    ! SETUP POINTERS TO LAGRANGIAN TRACER MEMORY AND TRACER INFORMATION
    CALL get_tracer_set(status, LGTRSTR, trlist_lg, ti_lg, ntrac_lg &
         , xt=pxt_a, xtte=pxtte_a, xtm1=pxtm1_a, xmem=pxtf_a)
    CALL tracer_halt(substr, status)

    til => trlist_lg
    DO
      IF (.NOT. ASSOCIATED(til)) EXIT
      IF (til%info%meta%cask_i(I_mix) == ON) THEN
        number_mix = number_mix + 1
      END IF
      til => til%next
    END DO

    conversion_lg2gp: IF (l_conv_lg2gp) THEN

       ! COPY META-INFORMATION
       CALL copy_tracer_set(status, LGTRSTR, LGGPTRSTR)
       CALL tracer_halt(substr, status)

       dims(1) = nproma
       dims(2) = nlev
       dims(3) = ngpblks

       ! only one time level, no standard, tracers shall not be externally
       ! initialized
       CALL setup_tracer_set(status, LGGPTRSTR, dims, 1, .FALSE., .FALSE.)
       CALL tracer_halt(substr, status)

       ! SETUP POINTERS TO MEMORY AND TRACER INFORMATION
       CALL get_tracer_set(status, LGGPTRSTR, trlist_lggp, ti_lggp &
            , ntrac_lggp, xt=pxt_lggp)
       CALL tracer_halt(substr, status)

    END IF conversion_lg2gp

  END SUBROUTINE setup_tracer_set_lg
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE setup_tracer_set_om
!!#D mpiom +
    ! SETUP MEMORY FOR TRACERS IN OCEAN

    USE messy_mpiom_mem_e5, ONLY: IE,JE,KE

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'setup_tracer_set_om'
    INTEGER                        :: status
    INTEGER, DIMENSION(3)          :: dims

    dims(1) = IE
    dims(2) = JE
    dims(3) = KE

    ! 4 time levels: xt, xtte, xtm1, xtf (the latter is 'extended')
    ! first 3 time levels are 'standard' (t, tendency, t-1) -> .TRUE.
    ! tracers can be initialiszed -> .TRUE.
    CALL setup_tracer_set(status, OMTRSTR, dims, 1, .TRUE., .TRUE.)
    CALL tracer_halt(substr, status)

    ! SETUP POINTERS TO LAGRANGIAN TRACER MEMORY AND TRACER INFORMATION
    CALL get_tracer_set(status, OMTRSTR, trlist_om, ti_om, ntrac_om &
         , xt=pxt_om )
    CALL tracer_halt(substr, status)
!!#D mpiom -
  END SUBROUTINE setup_tracer_set_om
  ! -------------------------------------------------------------------
#endif

! =====================================================================
#ifdef ICON
  ! -------------------------------------------------------------------
  SUBROUTINE map_tracer_meta_i2m(flag)

    USE messy_main_bmluse_bi,     ONLY: p_nh_state_lists, t_list_element &
                                      , t_var_metadata, t_var_metadata_dynamic &
                                      , get_var_name, n_dom
    USE messy_main_tools,         ONLY: int2str
    USE messy_main_blather_bi,    ONLY: error_bi, warning_bi, info_bi
    USE messy_main_data_bi,       ONLY: current_timelev
    USE messy_main_channel_mem,   ONLY: dom_current

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)           :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'map_tracer_meta_i2m'
    INTEGER :: jg
    INTEGER :: jt, it
    TYPE(t_list_element),             POINTER :: element
    TYPE (t_var_metadata),            POINTER :: info
    TYPE (t_var_metadata_dynamic),    POINTER :: info_dyn
    INTEGER :: status
    CHARACTER(LEN=3) :: tstr = '   '
    INTEGER :: zid
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER :: id_H2O

    SELECT CASE(flag)

    CASE(1) ! -----------------------------------------------------------------

       ! create empty dummy entries for ICON tracers (all patches)

       patch_loop1: DO jg=1, n_dom

          tracer_loop: DO jt=1, l_ntrac_icon(jg)

             CALL int2str(tstr, jt)

             CALL new_tracer(status, L_GPTRSTR(jg)    &
               , 'dummy'//tstr, 'icon'                &
               , unit = 'undefined'                   &
               )
             CALL tracer_halt(substr//'(1)', status)

          END DO tracer_loop

       END DO patch_loop1

    CASE(2) ! -----------------------------------------------------------------

       ! loop over icon tracer meta information (all patches) and
       ! update corresponding MESSy tracer meta information

       jg = dom_current

       ! tracer meta information different for diffrent time levels?
       it = current_timelev
       IF (it > 1) RETURN

       ! TRACER SET
       CALL get_tracer_set_id(status, L_GPTRSTR(jg), zid)
       CALL tracer_halt(substr//'(2)', status)

       element => &
            p_nh_state_lists(jg)%prog_list(it)%p%first_list_element
            ! t_nh_state_lists ! t_var_list |
            ! -> TYPE(t_var_list_intrinsic, POINTER  %p

       jt = 0
       element_loop: DO ! loop over elements in linked list
          element => element%next_list_element
          IF (.NOT.ASSOCIATED(element)) EXIT

          info     => element%field%info
          info_dyn => element%field%info_dyn

          IF (.NOT. info_dyn%tracer%lis_tracer) CYCLE

          jt = jt + 1
          CALL int2str(tstr, jt)

          IF (jt > l_ntrac_icon(jg)) THEN
             CALL error_bi('NUMBER OF TRACERS MISMATCH', substr//'(2)')
          ENDIF

          CALL info_bi( &
               ' mapping '//TRIM(get_var_name(element%field))//&
               &' meta-info; index = '//tstr &
               , substr//'(2)')

          ! SEARCH FOR EXISITNG TRACER IN LIST
          ti => trset(zid)%tilist
          DO
             IF (.NOT. ASSOCIATED(ti)) EXIT
             IF (TRIM(ti%info%ident%fullname) == 'dummy'//tstr) EXIT ! FOUND
             ti => ti%next
          END DO

          IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
             CALL error_bi('DUMMY TRACER NOT FOUND', substr//'(2)')
          END IF

          ti%info%ident%basename = TRIM(get_var_name(element%field))
          ti%info%ident%subname  = ''
          ti%info%ident%fullname = TRIM(ti%info%ident%basename)
          ti%info%ident%longname = TRIM(info%cf%long_name)
          ti%info%ident%standardname = TRIM(info%cf%standard_name)
          ti%info%ident%unit     = TRIM(info%cf%units)
          ti%info%ident%submodel = 'icon'
          !ti%info%ident%medium   = ! AIR
          !ti%info%ident%quantity = ! AMOUNTFRACTION
          !ti%info%ident%type     = ! SINGLE
          !
          ! ti%info%meta%cask_i
          ! ti%info%meta%cask_r
          ! ti%info%meta%cask_s

          ! avoid initialisation of ICON (q*) tracers as MESSy TRACERs
          ! with 0.0_dp ...
          ti%info%meta%linit = .TRUE. ! already initialsed

          ! set tracer properties of ICON hydrolocical cycle tracers
          ! to those of species 'H2O' ... see mo_nonhydro_state.f90, add_ref
          !
          id_H2O = get_chemprop_index('H2O')
          IF (id_H2O > 0) THEN
             SELECT CASE(TRIM(ti%info%ident%basename))
             CASE('qv','qc','qi','qr','qs','qg','qh')
                ti%info%meta%cask_i(1:MAX_CASK_I_CHEMPROP) = &
                     chemprop(id_H2O)%cask_i(:)
                ti%info%meta%cask_r(1:MAX_CASK_R_CHEMPROP) = &
                     chemprop(id_H2O)%cask_r(:)
                ti%info%meta%cask_s(1:MAX_CASK_S_CHEMPROP) = &
                     chemprop(id_H2O)%cask_s(:)
             ! for 2 moment micropysics scheme number concentrations
             CASE('qnc','qni','qnr','qns','qng','qnh','qninuc' &
                   ,'ninact','nccn','ninpot')
                ! set chemical tracer properties to values of water
                ti%info%meta%cask_i(1:MAX_CASK_I_CHEMPROP) = &
                     chemprop(id_H2O)%cask_i(:)
                ti%info%meta%cask_r(1:MAX_CASK_R_CHEMPROP) = &
                     chemprop(id_H2O)%cask_r(:)
                ti%info%meta%cask_s(1:MAX_CASK_S_CHEMPROP) = &
                     chemprop(id_H2O)%cask_s(:)
                ! overwrite molar mass
                ti%info%meta%cask_r(R_molarmass) = 1._dp
             CASE DEFAULT
                ! do nothing
             END SELECT
          END IF

       END DO element_loop

       IF (jt < l_ntrac_icon(jg)) THEN
          CALL warning_bi('NUMBER OF ICON TRACERS DECREASED', substr//'(2)')
       END IF

    END SELECT

  END SUBROUTINE map_tracer_meta_i2m
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE map_tracer_meta_m2i

    ! ICON
    USE mo_add_tracer_ref,       ONLY: add_tracer_ref
    USE mo_cf_convention,        ONLY: t_cf_var
    USE mo_grib2,                ONLY: t_grib2_var, grib2_var
    USE mo_cdi,                  ONLY: DATATYPE_FLT32, DATATYPE_FLT64
    USE mo_io_config,            ONLY: lnetcdf_flt64_output
    USE mo_cdi,                  ONLY: GRID_UNSTRUCTURED, DATATYPE_PACK16
    USE mo_cdi_constants,        ONLY: GRID_CELL
    USE mo_parallel_config,      ONLY: nproma
    USE mo_impl_constants,       ONLY: TLEV_NNOW_RCF
    USE mo_tracer_metadata,      ONLY: create_tracer_metadata
    USE mo_var_list,             ONLY: get_timelevel_string
    !
    ! ... to avoid parameter lists ...
    USE messy_main_bmluse_bi,    ONLY: p_nh_state_lists, p_nh_state &
                                     , p_patch &
                                     , advection_config

    ! MESSy
    USE messy_main_blather_bi,   ONLY: error_bi, warning_bi, info_bi
    USE messy_main_data_bi,      ONLY: current_timelev
    USE messy_main_channel_mem,  ONLY: dom_current
    USE messy_main_tools,        ONLY: int2str

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'map_tracer_meta_m2i'
    INTEGER           :: status
    INTEGER           :: jg
    INTEGER           :: it
    INTEGER           :: jt
    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc
    INTEGER           :: datatype_flt !< floating point accuracy in NetCDF output
    INTEGER           :: ibits        !< "entropy" of horizontal slice
    CHARACTER(LEN=2)  :: pstr = '  '
    CHARACTER(LEN=2)  :: tstr = '  '
    CHARACTER(LEN=3)  :: trstr = '   '
    !
    INTEGER           :: dummy_idx, nlev, nblks_c
    INTEGER           :: shape3d_c(3)
    CHARACTER(len=4)  :: suffix
    !
    ! ICON TRACER PROCESS SWITCHES AND TRACER PROPERTIES
    ! (see src/shared/mo_tracer_metadata.f90)
    INTEGER :: ihadv_tracer
    INTEGER :: ivadv_tracer
    LOGICAL :: lturb_tracer
    LOGICAL :: lconv_tracer
    !INTEGER :: ised_tracer
    !LOGICAL :: ldep_tracer
    !INTEGER :: iwash_tracer
    !REAL(DP) :: mol_weight

    jg = dom_current
    it = current_timelev

    CALL int2str(pstr, jg, ' ')
    CALL int2str(tstr, it, ' ')

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF
    ibits = DATATYPE_PACK16

    nlev      = p_patch(jg)%nlev
    nblks_c   = p_patch(jg)%nblks_c
    shape3d_c = (/ nproma, nlev, nblks_c /)
    suffix    = get_timelevel_string(it)

    CALL get_tracer_set(status, L_GPTRSTR(jg), ti = ti_gp, ntrac = ntrac_gp)
    CALL tracer_halt(substr, status)

    DO jt = ntrac_icon+1, ntrac_gp, 1

       CALL int2str(trstr, jt, ' ')

       CALL info_bi('  '//TRIM(ti_gp(jt)%tp%ident%fullname)//&
            &' (idt = '//TRIM(trstr)//'; '//&
            &' patch = '//pstr//'; timelev = '//tstr//')'   &
            , substr)

       ! MAP MESSy TRACER PROPERTIES TO ICON PROPERTIES

       ! see mo_advection_config.f90
       SELECT CASE (ti_gp(jt)%tp%meta%cask_i(I_advect))
       CASE(OFF)
          !
          ihadv_tracer = 0
          ivadv_tracer = 0
          !
       CASE(ON)
          !
          ! ICON defaults
          !ihadv_tracer = 2
          !ivadv_tracer = 3
          !
          ihadv_tracer = 22 ! 52
          ivadv_tracer = 3
          !
       CASE DEFAULT
          !
          ! 1000*ihadv + ivadv
          ihadv_tracer = ti_gp(jt)%tp%meta%cask_i(I_ADVECT) / 1000
          ivadv_tracer = MOD(ti_gp(jt)%tp%meta%cask_i(I_ADVECT), 1000)
          !
       END SELECT

       !lfeedback (from child to parent domain)

       lturb_tracer = ti_gp(jt)%tp%meta%cask_i(I_VDIFF)   == ON
       lconv_tracer = ti_gp(jt)%tp%meta%cask_i(I_CONVECT) == ON
       !ised_tracer
       !ldep_tracer
       !iwash_tracer
       !
       !solubility
       !rho
       !mol_weight = ti_gp(jt)%tp%meta%cask_r(R_molarmass)

       cf_desc    = t_cf_var(TRIM(ti_gp(jt)%tp%ident%fullname)    &
            , TRIM(ti_gp(jt)%tp%ident%unit)                       &
            , 'MESSy tracer', datatype_flt)
       grib2_desc = grib2_var(255,255,255,ibits,GRID_UNSTRUCTURED,GRID_CELL)

       IF (it == 1) THEN
       CALL add_tracer_ref( p_nh_state_lists(jg)%prog_list(it), 'tracer' &
            , TRIM(ti_gp(jt)%tp%ident%fullname)//suffix     &
            , dummy_idx                             &
            , p_nh_state(jg)%prog(it)%tracer_ptr(:) &
            , cf_desc, grib2_desc                   &
            , advection_config(p_patch(jg)%id)      &
            , jg          = p_patch(jg)%id          &
            , ldims       = shape3d_c               &
            , loutput     = .TRUE.                  &
            , lrestart    = .TRUE.                  &
            , tlev_source = TLEV_NNOW_RCF           &
!!$         , tlev_source = it           &
            , tracer_info = create_tracer_metadata( lis_tracer=.TRUE. &
            , lfeedback   = ( ti_gp(jt)%tp%meta%cask_i(I_FEEDBACK) == ON ) &
            , name = TRIM(ti_gp(jt)%tp%ident%fullname)//suffix &
            , ihadv_tracer = ihadv_tracer &
            , ivadv_tracer = ivadv_tracer &
            , lturb_tracer = lturb_tracer &
            , lconv_tracer = lconv_tracer &
!            , ised_tracer  = ised_tracer &
!            , ldep_tracer  = ldep_tracer &
!            , iwash_tracer = iwash_tracer &
            ) )
       ELSE
       CALL add_tracer_ref( p_nh_state_lists(jg)%prog_list(it), 'tracer' &
            , TRIM(ti_gp(jt)%tp%ident%fullname)//suffix     &
            , dummy_idx                             &
            , p_nh_state(jg)%prog(it)%tracer_ptr(:) &
            , cf_desc, grib2_desc                   &
            , advection_config(p_patch(jg)%id)      &
            , jg          = p_patch(jg)%id          &
            , ldims       = shape3d_c               &
            , loutput     = .TRUE.                  &
            , lrestart    = .TRUE.                  &
            , tlev_source = TLEV_NNOW_RCF           &
!!$         , tlev_source = it           &
            )
       END IF

       advection_config(p_patch(jg)%id)%itype_vlimit(jt) = 1 !2
       advection_config(p_patch(jg)%id)%itype_hlimit(jt) = 3

!qqq translate further MESSy meta information switches (cask_I;R;C)
!    (I_ADV, I_CONV, .. ...)
!    for create_tracer_metadata
! NOTE: _hydro is empty and _chem only relevant for ART

    END DO

  END SUBROUTINE map_tracer_meta_m2i
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE map_tracer_memory_icon(flag)

    ! ICON
    ! ... to avoid parameter lists ...
    USE mo_nonhydro_state,       ONLY: p_nh_state_lists, p_nh_state
    USE mo_var_list,             ONLY: find_list_element
    USE mo_linked_list,          ONLY: t_list_element
    USE mo_kind,                 ONLY: i8
    USE mo_nwp_phy_state,        ONLY: prm_nwp_tend
    USE messy_main_bmluse_bi,    ONLY: p_patch, nnew, nnow

    ! MESSy
    USE messy_main_channel_bi,      ONLY: n_dom
    USE messy_main_data_bi,         ONLY: current_timelev
    USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev, ngpblks
    USE messy_main_channel_mem,  ONLY: dom_current
    USE messy_main_blather_bi,   ONLY: info_bi, error_bi
    USE messy_main_tools,        ONLY: int2str
    USE messy_main_tendency_bi,  ONLY: mtend_get_handle &
                                     , mtend_register   &
                                     , mtend_request

    IMPLICIT NONE
    INTRINSIC :: INT, PRODUCT

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'map_tracer_memory_icon'
    INTEGER :: jg, jt, status, it
    TYPE(t_list_element), POINTER :: list_element
    INTEGER :: ndims
    INTEGER :: idims(5)
    CHARACTER(LEN=3)  :: ctrstr = '   '
    CHARACTER(LEN=3)  :: ttrstr = '   '
    INTEGER :: nconv_tracer, nturb_tracer
    INTEGER :: ic, jc, jb, nblks_c
    INTEGER :: conv_handle, turb_handle
    REAL(dp), POINTER, DIMENSION(:,:,:) :: tend_ptr => NULL()

    jg = dom_current
    ! current_timelev is the working timelev (nothing to do with running time)
    it = current_timelev
    nblks_c = p_patch(dom_current)%nblks_c

    CALL get_tracer_set(status, L_GPTRSTR(jg) &
         , trlist_gp, ti_gp, ntrac_gp         &
         , xt=pxt, xtte=pxtte, xtm1=pxtm1)
    CALL tracer_halt(substr, status)

    list_element => find_list_element(p_nh_state_lists(jg)%prog_list(it) &
         , 'tracer')

    ! TRACERS FOR CONVECTIVE TRANSPORT
    ! Note: art_config(jg)%nconv_tracer and %nturb_tracer are
    !       not yet up to date; they are counted later ...
    !       Therefore, we do our own counting here.
    conv_handle=mtend_get_handle('conv')

    nconv_tracer = 0
    DO jt = ntrac_icon+1, ntrac_gp, 1
       IF (ti_gp(jt)%tp%meta%cask_i(I_CONVECT) == ON ) &
            nconv_tracer = nconv_tracer + 1
    END DO
    !
    CALL int2str(ctrstr, nconv_tracer, ' ')
    CALL info_bi('  number of tracers for convective transport: '//ctrstr &
         , substr)

    ! TRACERS FOR TURBULENT TRANSPORT
    nturb_tracer = 0
    DO jt = ntrac_icon+1, ntrac_gp, 1
       IF (ti_gp(jt)%tp%meta%cask_i(I_VDIFF) == ON ) &
            nturb_tracer = nturb_tracer + 1
    END DO
    !
    CALL int2str(ttrstr, nturb_tracer, ' ')
    CALL info_bi('  number of tracers for turbulent transport : '//ttrstr &
         , substr)

    SELECT CASE(flag)

    CASE(1) ! ---------------------------------------------------------
       ! deallocate ICON tracer pointers
       DEALLOCATE(list_element%field%r_ptr)
       NULLIFY(list_element%field%r_ptr)
       ! avoid deallocation in finalising phase; see mo_var_list.f90
       list_element%field%info%allocated = .FALSE.
       ! maybe not required ...
       NULLIFY(list_element%field%s_ptr)
       NULLIFY(list_element%field%i_ptr)
       NULLIFY(list_element%field%l_ptr)
       ! correct memory size counter to avoid error in finalising phase
       ! see mo_var_list.f90
       ndims = list_element%field%info%ndims
       idims(1:ndims)    = list_element%field%info%used_dimensions(1:ndims)
       idims((ndims+1):) = 1
       p_nh_state_lists(jg)%prog_list(it)%p%memory_used =      &
            p_nh_state_lists(jg)%prog_list(it)%p%memory_used - &
            INT(list_element%field%var_base_size, i8) *        &
            INT(PRODUCT(idims(1:5)),i8)

       ! point with ICON pointers to the corresponding MESSy TRACER memory
!qqq
       ! (see setup_tracer_set_gp for ICON ....)
       !QQQ do we need nnew/ nnow or nnew_rcf /nnow_rcf here?
       IF (it == nnew(jg)) THEN
          list_element%field%r_ptr => pxt(:,:,:,:,:)
       ELSE IF (it == nnow(jg)) THEN
          list_element%field%r_ptr => pxtm1(:,:,:,:,:)
       ELSE
          CALL error_bi('ERROR in tracer time levels init_memory(11)', substr)
       END IF

       p_nh_state(jg)%prog(it)%tracer => list_element%field%r_ptr(:,:,:,:,1)

       ! CONVECTIVE TRACERS
       ALLOCATE(p_nh_state(jg)%prog(it)%conv_tracer(nblks_c, nconv_tracer))
       !
       ic = 0
       conv_handle=mtend_get_handle('conv')
       DO jt = ntrac_icon+1, ntrac_gp, 1
          IF (ti_gp(jt)%tp%meta%cask_i(I_CONVECT) == ON ) THEN
             ic = ic + 1
             CALL mtend_register(conv_handle, jt, dom_id=jg)
             DO jb = 1, nblks_c
                p_nh_state(jg)%prog(it)%conv_tracer(jb,ic)%idx_tracer = jt
                ! nproma,nlev,ngpblks,ntrac_gp,1
                p_nh_state(jg)%prog(it)%conv_tracer(jb,ic)%ptr => &
!                     list_element%field%r_ptr(:,:,jb,jt,1)
! convection always accesses nnew, therefore always link nnew to avoid swapping
! is this correct / possible or causes this errors in other places?
                     pxt(:,:,jb,jt,1)
             END DO
          END IF
       END DO

       ! TURBULENT TRACERS
       ALLOCATE(p_nh_state(jg)%prog(it)%turb_tracer(nblks_c, nturb_tracer))
       ic = 0
       turb_handle=mtend_get_handle('turbdiff')

       DO jt = ntrac_icon+1, ntrac_gp, 1
          IF (ti_gp(jt)%tp%meta%cask_i(I_VDIFF) == ON ) THEN
             ic = ic + 1
             CALL mtend_register(turb_handle, jt, dom_id=jg)
             DO jb = 1, nblks_c
                p_nh_state(jg)%prog(it)%turb_tracer(jb,ic)%idx_tracer = jt
                ! nproma,nlev,ngpblks,ntrac_gp,1
                p_nh_state(jg)%prog(it)%turb_tracer(jb,ic)%ptr => &
!                     list_element%field%r_ptr(:,:,jb,jt,1)
                     pxt(:,:,jb,jt,1)
             END DO
          END IF
       END DO

       ! ALLOCATE air mass help fields
       IF (jg == 1) THEN
          ALLOCATE(grmass_conv(n_dom))
          ALLOCATE(grmassdry_conv(n_dom))
       END IF
       ALLOCATE(grmass_conv(jg)%ptr(nproma,nlev,nblks_c))
       ALLOCATE(grmassdry_conv(jg)%ptr(nproma,nlev,nblks_c))
       grmass_conv(jg)%ptr = 0._dp
       grmassdry_conv(jg)%ptr = 0._dp

    CASE(2) ! ---------------------------------------------------------

       ! NOTES:
       ! (1) prm_nwp_tend is not yet allocated at entry point 1;
       !     therefore we need a second call here to allocate/associate
       !     the memory for the convective tracer tendencies.
       ! (2) The convective tracer tendencies are NOT purely diagnostic.
       !     The call sequence is (2.4.0):
       !     - mo_nwp_conv_interface::nwp_convection
       !       CALL cumastrn(..., ktrac(IN), pcen(IN), ptenc(INOUT))
       !     - mo_cumastr::cumastrn:
       !       CALL cuctracer(...,ktrac(IN), ..., pcen(IN), ptenc(INOUT))
       !     - mo_cuflxtends::cuctracer: ...
       !     Thus, tracers are NOT directly modified by convection,
       !     i.e., the tracer tendencies need to be added explicitely.
       !     This is done (FOR ART ONLY) in mo_util_phys::tracer_add_phytend.
       ! (3) Strategy for MESSy:
       !     + ALLOCATE the memory for convective tracer tendencies
       !       within the ICON structure, in order to avoid additional
       !       MESSy-hooks in mo_nwp_conv_interface (call of cumastrn).
       !     + Make sure that art_config - structure contains correct
       !       entries. This should automatically happen by using
       !       add_tracer_ref in map_tracer_meta_m2i above.
       !     + Free the memory for convective tracer tendencies during
       !       finalize!
       !     + Tracers are supposed to be in kg/kg(mois air?) for convection.
       !       Thus conversion routines before and after the convection need
       !       to be called (via messy_convec(1) and messy_convec(2))
       !       which call the conversions: mol/mol(dry) <-> kg/kg(moist).

       ! TENDENCIES
       ALLOCATE(prm_nwp_tend(jg)%conv_tracer_tend(nblks_c, nconv_tracer))
       !
       ic = 0
       DO jt = ntrac_icon+1, ntrac_gp, 1
          IF (ti_gp(jt)%tp%meta%cask_i(I_CONVECT) == ON ) THEN
             ic = ic + 1
             CALL mtend_request('conv',jt,tend_ptr, jg, lextern=.TRUE.)
             DO jb = 1, nblks_c
                prm_nwp_tend(jg)%conv_tracer_tend(jb, ic)%idx_tracer = jt
                ! Note:
                ! Pointing to pxtte, the MESSy tracer tendency would
                ! be overwritten with the convective tracer tendency.
                ! All other tendencies would be lost.
                ! Therefore, allocate the ICON internal memory,
                ! which also reduces required code changes in ICON itself.
                !ALLOCATE(&
                !  prm_nwp_tend(jg)%conv_tracer_tend(jb,ic)%ptr(nproma,nlev) )
                prm_nwp_tend(jg)%conv_tracer_tend(jb,ic)%ptr => tend_ptr(:,:,jb)
                prm_nwp_tend(jg)%conv_tracer_tend(jb, ic)%ptr(:,:) = 0.0_dp
                !
             END DO
             NULLIFY(tend_ptr)
          END IF
       END DO

       ! TENDENCIES
       ALLOCATE(prm_nwp_tend(jg)%turb_tracer_tend(nblks_c, nturb_tracer))
       ic = 0
       DO jt = ntrac_icon+1, ntrac_gp, 1
          IF (ti_gp(jt)%tp%meta%cask_i(I_VDIFF) == ON ) THEN
             ic = ic + 1
             CALL mtend_request('turbdiff',jt,tend_ptr, jg, lextern=.TRUE.)
             DO jb = 1, nblks_c
                prm_nwp_tend(jg)%turb_tracer_tend(jb, ic)%idx_tracer = jt
                ! Note:
                ! Pointing to pxtte, the MESSy tracer tendency would
                ! be overwritten with the convective tracer tendency.
                ! All other tendencies would be lost.
                ! Therefore, allocate the ICON internal memory,
                ! which also reduces required code changes in ICON itself.
                prm_nwp_tend(jg)%turb_tracer_tend(jb,ic)%ptr => tend_ptr(:,:,jb)
                prm_nwp_tend(jg)%turb_tracer_tend(jb,ic)%ptr(:,:) = 0.0_dp
                !
             END DO
             NULLIFY(tend_ptr)
          END IF
       END DO
    END SELECT

  END SUBROUTINE map_tracer_memory_icon
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_turbdiff(flag, noffset, turbptr)

    USE messy_main_tendency_bi,     ONLY: mtend_get_handle   &
                                        , mtend_add_l
    USE messy_main_grid_def_mem_bi, ONLY: jrow, nlev, nproma &
                                        , startidx, endidx
    USE messy_main_grid_def_bi,     ONLY: grmass, grmassdry   &
                                       , update_airmass
    USE messy_main_data_bi,         ONLY: current_timelev
    USE messy_main_channel_mem,     ONLY: dom_current
    USE messy_main_bmluse_bi,       ONLY: modvar, atm_phy_nwp_config &
                                        , prm_nwp_tend

    IMPLICIT NONE

    ! I/O
    INTEGER,      INTENT(IN)                   :: flag
    INTEGER,      INTENT(IN)                   :: noffset
    TYPE(modvar), INTENT(INOUT), DIMENSION(:)  :: turbptr

    ! LOCAL
    INTEGER                             :: jg, jt, jk, jc, ic
    REAL(DP), DIMENSION(:,:,:), POINTER :: ptr => NULL()
    REAL(DP)                            :: pdtime
    REAL(DP), DIMENSION(:,:),   POINTER :: tend => NULL()
    INTEGER                             :: turb_handle

    jg = dom_current
    pdtime = atm_phy_nwp_config(jg)% dt_fastphy

    NULLIFY(ptr)
    ptr => pxt(:,:,jrow,:,1)
    SELECT CASE(flag)

    CASE (1)
       !
       ! Note: MESSy TRACERs are in mol/mol(dry) and need to be converted
       !       into kg/kg(moist) for the call to convectionturbulence.
       !
       ! Note: grmass and grmassdry need to be updated here!
       CALL update_airmass(jrow)
       !
      ! QQQQtest if same as in convection
       !
       ic = 0
       DO jt = ntrac_icon+1, ntrac_gp, 1
          IF (ti_gp(jt)%tp%meta%cask_i(I_VDIFF) == ON ) THEN
             ic = ic + 1
             !
             ! Note: The scaling with constant molar mass is omitted here,
             !       becasue the backward conversion (case 2) would
             !       require the inverse scaling ...
             !
             !       kg(tracer)   mol(tracer)    M_tracer    grmass_dry
             !       ---------- = ----------- (* --------) * ----------
             !       kg(moist)    mol(dry)       M_dryair    grmass_moist
             !
             DO jk = 1, nlev
                DO jc = startidx, endidx
                   ptr(jc,jk,jt) = ptr(jc,jk,jt) &
                        * grmassdry(jc,jk,jrow)/grmass(jc,jk,jrow)
                END DO
             END DO

             turbptr(noffset+ic)%av => pxt(:,:,jrow,jt,1)
             turbptr(noffset+ic)%at => &
                  prm_nwp_tend(jg)%turb_tracer_tend(jrow,ic)%ptr(:,:)
             turbptr(noffset+ic)%at =  0._dp
             turbptr(noffset+ic)%fc = .FALSE.
          END IF
       END DO

    CASE (2)
       ! FIRST apply tendencies : TODO add to the xtte instead ?
       ! SECOND CONVERT TRACER BACK to mol/mol(dry)
       ic = 0
       turb_handle=mtend_get_handle('turbdiff')
       ALLOCATE(tend(nproma,nlev))
       DO jt = ntrac_icon+1, ntrac_gp, 1
          IF (ti_gp(jt)%tp%meta%cask_i(I_VDIFF) == ON ) THEN
             ic = ic + 1
             DO jk = 1,nlev
                DO jc = startidx, endidx
!!$                   ! 1) This would be the same implementation as in the ICON
!!$                   turbptr(noffset+ic)%av(jc,jk) = &
!!$                        (MAX(0._dp,turbptr(noffset+ic)%av(jc,jk)  &
!!$                        + pdtime * turbptr(noffset+ic)%at(jc,jk))
!!$
                   tend(jc,jk) = turbptr(noffset+ic)%at(jc,jk) &
                        * grmass(jc,jk,jrow)/ grmassdry(jc,jk,jrow)
                   ! 2)
                   ptr(jc,jk,jt) = ptr(jc,jk,jt)  &
                        * grmass(jc,jk,jrow)/grmassdry(jc,jk,jrow)
                   !
                END DO
             END DO
             CALL mtend_add_l(turb_handle,jt, px=tend)
          END IF
       END DO
       DEALLOCATE(tend); NULLIFY(tend)
    END SELECT


  END SUBROUTINE main_tracer_turbdiff
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_convec(flag)

    USE messy_main_tendency_bi,     ONLY: mtend_get_handle   &
                                        , mtend_add_l
    USE messy_main_grid_def_mem_bi, ONLY: jrow, nlev, nproma &
                                        , startidx, endidx
    USE messy_main_grid_def_bi,     ONLY: grmass, grmassdry
    USE messy_main_data_bi,         ONLY: current_timelev
    USE messy_main_channel_mem,     ONLY: dom_current


    ! ICON
    USE messy_main_bmluse_bi, ONLY: prm_nwp_tend, atm_phy_nwp_config

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    INTEGER                             :: jg, jt, jk, jc, ic
    REAL(DP), DIMENSION(:,:,:), POINTER :: ptr => NULL()
    REAL(DP), DIMENSION(:,:),   POINTER :: zm  => NULL()
    REAL(DP), DIMENSION(:,:),   POINTER :: zmd => NULL()
    REAL(DP)                            :: zte
    REAL(DP)                            :: pdtime
    REAL(DP), DIMENSION(:,:),   POINTER :: tend => NULL()
    INTEGER                             :: conv_handle

    jg = dom_current
    pdtime = atm_phy_nwp_config(jg)% dt_fastphy

    NULLIFY(ptr)
    ptr => pxt(:,:,jrow,:,1)

    SELECT CASE(flag)
       !
    CASE(1)
       !
       ! Note: MESSy TRACERs are in mol/mol(dry) and need to be converted
       !       into kg/kg(moist) for the call to convection.
       !
       ! Note: grmass and grmassdry need to be updated here!
       ! 2020-10-21:
       ! most probably this call to update_airmassis not required (anymore),
       ! as it was already called in messy_turbdiff.... comparing the grmasses
       ! leads to identical results, might become necessary later on, if more
       ! processes are coupled
       !
       ! CALL update_airmass(jrow)
       !
       grmass_conv(jg)%ptr(:,:,jrow)    = grmass(:,:,jrow)
       grmassdry_conv(jg)%ptr(:,:,jrow) = grmassdry(:,:,jrow)
       !
       ic = 0
       DO jt = ntrac_icon+1, ntrac_gp, 1
          IF (ti_gp(jt)%tp%meta%cask_i(I_CONVECT) == ON ) THEN
             ic = ic + 1
             !
             ! Note: The scaling with constant molar mass is omitted here,
             !       becasue the backward conversion (case 2) would
             !       require the inverse scaling ...
             !
             !       kg(tracer)   mol(tracer)    M_tracer    grmass_dry
             !       ---------- = ----------- (* --------) * ----------
             !       kg(moist)    mol(dry)       M_dryair    grmass_moist
             !
             DO jk = 1, nlev
                DO jc = startidx, endidx
                   ptr(jc,jk,jt) = ptr(jc,jk,jt) &
                        * grmassdry(jc,jk,jrow)/grmass(jc,jk,jrow)
                END DO
             END DO
             ! reset tracer tendency before new convection calculation
             prm_nwp_tend(jg)%conv_tracer_tend(jrow,ic)%ptr(:,:) = 0._dp
          END IF
       END DO
       !
    CASE(2)
       !
       ! Note: The TRACERs do not have changed by convection, only
       !       the tendencies have been calculated. These are, however,
       !       in kg/kg(moist)/s. For (later) use of MESSy TENDENCY,
       !       they need to be converted into mol/mol(dry)/s.
       !       Since the (moist) airmass has been changed by convection,
       !       for numerical consistency the tendencies need to be recalculated
       !       from the total change ...
       !       For this conversion, we need to save the "old" airmasses
       !       to convert back (into mol/mol(dry)) the unchanged TRACERs.
       !
       ic = 0
       DO jt = ntrac_icon+1, ntrac_gp, 1
          IF (ti_gp(jt)%tp%meta%cask_i(I_CONVECT) == ON ) THEN
             ic = ic + 1

             DO jk = 1, nlev
                DO jc = startidx, endidx
                   !
                   ! convert back unchanged TRACER
                   ptr(jc,jk,jt) = ptr(jc,jk,jt)  &
                        * grmass(jc,jk,jrow)/grmassdry(jc,jk,jrow)
                   !
                END DO
             END DO
          END IF
       END DO
       !
    CASE (3)

       ic = 0
       conv_handle=mtend_get_handle('conv')
       ALLOCATE(tend(nproma,nlev))
       ! (3) Rescale CONVECTIVE TENDENCY ACCORDING TO DRYMASS CHANGES
       !     Note: mass is only conserved, if we apply the tendency here
       !           directly to the current (nnew) tracer field (xt).

       ! The memory of TENDENCY is used for conv_tracer_tend to reduce the
       ! memory food print. However, to also use the accounting system of
       ! TENDENCY correctly, we have to set the tendency by call of
       ! mtend_add_l/g. Therefore, here conv_tracer_tend is copied to a local
       ! value and set to zero. The local field is scaled according to the
       ! grid mass changes and finally, the new scaled tendency is assign in
       ! the TENDENCY conform way by call of mtend_add
       DO jt = ntrac_icon+1,ntrac_gp

          IF (ti_gp(jt)%tp%meta%cask_i(I_CONVECT) == ON ) THEN
             ic = ic + 1
             tend(:,:) = 0._dp
             DO jc = startidx, endidx
                DO jk = 1, nlev
                   IF (grmassdry_conv(jg)%ptr(jc,jk,jrow) > 0._dp) THEN

                      tend(jc,jk) = &
                          pdtime*prm_nwp_tend(jg)%conv_tracer_tend(jrow,ic)%ptr(jc,jk)&
                       * grmass_conv(jg)%ptr(jc,jk,jrow) &
                       / grmassdry_conv(jg)%ptr(jc,jk,jrow)

                   xt(_RI_XYZN_(jc,jrow,jk,jt)) = xt(_RI_XYZN_(jc,jrow,jk,jt)) &

                        + tend(jc,jk)

!!$                   xt(_RI_XYZN_(jc,jrow,jk,jt)) = xt(_RI_XYZN_(jc,jrow,jk,jt)) &
!!$
!!$                        + pdtime*prm_nwp_tend(jg)%conv_tracer_tend(jrow,ic)%ptr(jc,jk)&
!!$                       * grmass_conv(jg)%ptr(jc,jk,jrow)/grmassdry_conv(jg)%ptr(jc,jk,jrow)
                   END IF
                END DO
!!$               WHERE (grmassdry_conv(jg)%ptr(jc,:,jrow) > 0._dp)
!!$                   xt(_RI_XYZN_(jc,jrow,:,jt)) = xt(_RI_XYZN_(jc,jrow,:,jt)) &
!!$
!!$                        + pdtime*prm_nwp_tend(jg)%conv_tracer_tend(jrow,ic)%ptr(jc,:)&
!!$                        * grmass_conv(jg)%ptr(jc,:,jrow)/grmassdry_conv(jg)%ptr(jc,:,jrow)
!!$                END WHERE
             END DO
             tend(:,:) = tend(:,:) / pdtime
             CALL mtend_add_l(conv_handle, jt, px=tend, l_add=.FALSE.)
          END IF
       END DO
       DEALLOCATE(tend); NULLIFY(tend)

    END SELECT

  END SUBROUTINE main_tracer_convec
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE icon_tracer_integration

    USE messy_main_timer, ONLY: delta_time

    IMPLICIT NONE

    ! add MESSy tracer tendencies of MESSy submodels (for current patch)
    xt(:,:,:,:) = xt(:,:,:,:) + xtte(:,:,:,:)*delta_time
!    xtm1(:,:,:,:) = xtm1(:,:,:,:) + xtte(:,:,:,:)*delta_time

  END SUBROUTINE icon_tracer_integration
  ! -------------------------------------------------------------------
#endif
! =====================================================================

  ! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_fconv_loc(direction, callstr, TRSETSTR)

    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_grid_def_mem_bi,  ONLY: jrow, kproma
    USE messy_main_timer,            ONLY: time_step_len
    USE messy_main_tracer_family_bi, ONLY: tracfamily_1_f2t, tracfamily_1_t2f &
                                         , tracfamily_2_sum, tracfamily_2_rsc

    IMPLICIT NONE

    ! I/O
    CHARACTER(len=3), INTENT(in) :: direction
    CHARACTER(len=*), INTENT(in) :: callstr
    CHARACTER(len=*), INTENT(in) :: TRSETSTR

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_fconv_loc'
    INTEGER                     :: status

    IF (.NOT. l_family) RETURN

    SELECT CASE (direction)
       !
    CASE ('f2t','F2t','f2T','F2T')
       !
       CALL tracfamily_1_f2t(status, callstr, p_pe, TRSETSTR, &
            time_step_len, jrow, kproma)
       CALL tracer_halt(substr, status)
       !
    CASE ('t2f','T2f','t2F','T2F')
       !
       CALL tracfamily_1_t2f(status, callstr, p_pe, TRSETSTR, &
            time_step_len, jrow, kproma)
       CALL tracer_halt(substr, status)
       !
    CASE ('sum','SUM')
       !
       CALL tracfamily_2_sum(TRSETSTR, jrow)
       !
    CASE ('rsc','RSC')
       !
       CALL tracfamily_2_rsc(TRSETSTR, time_step_len, jrow)
       !
    CASE default
       !
       status = 2010 ! UNKNOWN CONVERSION FLAG
       CALL tracer_halt(substr, status)
       !
    END SELECT

  END SUBROUTINE main_tracer_fconv_loc
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_fconv_glb(direction, callstr, TRSETSTR)

    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, npromz, ngpblks
    USE messy_main_timer,            ONLY: time_step_len
    USE messy_main_tracer_family_bi, ONLY: tracfamily_1_f2t, tracfamily_1_t2f &
                                         , tracfamily_2_sum, tracfamily_2_rsc

    IMPLICIT NONE

    ! I/O
    CHARACTER(len=3), INTENT(in) :: direction
    CHARACTER(len=*), INTENT(in) :: callstr
    CHARACTER(len=*), INTENT(in) :: TRSETSTR

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_fconv_glb'
    INTEGER :: status
    INTEGER :: jjrow, jp

    IF (.NOT. l_family) RETURN

    SELECT CASE (direction)
       !
    CASE ('f2t','F2t','f2T','F2T')
       !
       DO jjrow = 1, ngpblks
#ifndef CESM1
          IF (jjrow == ngpblks) THEN
             jp = npromz
          ELSE
             jp = nproma
          END IF
#else
          jp = npromz(jjrow)
#endif
          CALL tracfamily_1_f2t(status, callstr, p_pe, TRSETSTR, &
               time_step_len, jjrow, jp)
          CALL tracer_halt(substr, status)
       END DO
       !
    CASE ('t2f','T2f','t2F','T2F')
       !
       DO jjrow = 1, ngpblks
#ifndef CESM1
          IF (jjrow == ngpblks) THEN
             jp = npromz
          ELSE
             jp = nproma
          END IF
#else
          jp = npromz(jjrow)
#endif
          CALL tracfamily_1_t2f(status, callstr, p_pe, TRSETSTR, &
               time_step_len, jjrow, jp)
          CALL tracer_halt(substr, status)
       END DO
       !
    CASE ('sum','SUM')
       !
       DO jjrow = 1, ngpblks
          CALL tracfamily_2_sum(TRSETSTR, jjrow)
       END DO
       !
    CASE ('rsc','RSC')
       !
       DO jjrow = 1, ngpblks
          CALL tracfamily_2_rsc(TRSETSTR, time_step_len, jjrow)
       END DO
       !
    CASE default
       !
       status = 2010 ! UNKNOWN CONVERSION FLAG
       CALL tracer_halt(substr, status)
       !
    END SELECT

  END SUBROUTINE main_tracer_fconv_glb
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_write_output(flag)

    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, npromz, ngpblks
    USE messy_main_timer,            ONLY: time_step_len
    USE messy_main_tracer_family_bi, ONLY: tracfamily_1_t2f

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_write_output'
    INTEGER :: status
    INTEGER :: jjrow, jp

!!#D attila +
    IF (l_conv_lg2gp) CALL tracer_out_conv_lg2gp(flag)
!!#D attila -
!!#D clams +
    IF (l_conv_cl2gp) CALL tracer_out_conv_cl2gp(flag)
!!#D clams -

    IF ( (l_family) .AND. (flag==1) ) THEN

       ! UPDATE FAMILIES TO BE CONSISTENT WITH TRACERS (FOR OUTPUT)
       DO jjrow = 1, ngpblks
#ifndef CESM1
          IF (jjrow == ngpblks) THEN
             jp = npromz
          ELSE
             jp = nproma
          END IF
#else
          jp = npromz(jjrow)
#endif
          CALL tracfamily_1_t2f(status, substr, p_pe, GPTRSTR, &
               time_step_len, jjrow, jp, .FALSE.)
          CALL tracer_halt(substr, status)
       END DO

    END IF

  END SUBROUTINE main_tracer_write_output
  ! ----------------------------------------------------------------------

!!#D attila +
  ! ----------------------------------------------------------------------
  SUBROUTINE tracer_out_conv_lg2gp(flag)

    USE messy_main_blather_bi,    ONLY: error_bi
#if defined(ECHAM5)
    USE messy_attila_tools_e5,    ONLY: lg2gp_e5
#endif

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tracer_out_conv_lg2gp'

    IF (ntrac_lggp /= ntrac_lg) &
         CALL error_bi('number of tracers mismatch',substr)

    SELECT CASE(flag)
    CASE(1)
       !
#if defined(ECHAM5)
       CALL lg2gp_e5(qxt_a, xt_lggp          &
            , i_conv_lg2gp_mode              &
            , fill_value = r_conv_lg2gp_fill &
            , lmcons = l_conv_lg2gp_mc       &
            )
#endif
       !
    CASE(2)
       !
!!$       xt_lggp(:,:,:,:) = 0.0_DP
       !
    CASE DEFAULT
       !
       CALL error_bi('UNKNOWN FLAG !', substr)

    END SELECT

  END SUBROUTINE tracer_out_conv_lg2gp
  ! ----------------------------------------------------------------------
!!#D attila -

!!#D clams +
  ! ----------------------------------------------------------------------
  SUBROUTINE tracer_out_conv_cl2gp(flag)

    USE messy_main_blather_bi,    ONLY: error_bi
#if defined(ECHAM5)
    USE messy_clams_tools_e5,     ONLY: cl2gp_e5
#endif
!!$    USE messy_main_constants_mem, ONLY: FLAGGED_BAD

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tracer_out_conv_cl2gp'

    IF (ntrac_clgp /= ntrac_cl) &
         CALL error_bi('number of tracers mismatch',substr)

    SELECT CASE(flag)
    CASE(1)
       !
#if defined(ECHAM5)
       CALL cl2gp_e5(qxt_c, xt_clgp          &
            , i_conv_cl2gp_mode              &
            , fill_value = r_conv_cl2gp_fill &
            )
#endif
       !
    CASE(2)
       !
!!$       xt_lggp(:,:,:,:) = 0.0_DP
       !
    CASE DEFAULT
       !
       CALL error_bi('UNKNOWN FLAG !', substr)

    END SELECT

  END SUBROUTINE tracer_out_conv_cl2gp
  ! ----------------------------------------------------------------------
!!#D clams -

  ! -------------------------------------------------------------------
  SUBROUTINE init_pvar(vname, dat)

    USE messy_main_timer,          ONLY: lresume
    USE messy_main_blather_bi,     ONLY: warning_bi, error_bi, info_bi
    USE messy_main_constants_mem,  ONLY: M_air, M_H2O
    USE messy_main_tools,          ONLY: int2str
    USE messy_main_data_bi,        ONLY: bmstr => modstr

#ifdef ECHAM5
   USE messy_main_data_bi,         ONLY: eps, q, qm1
#endif
#ifdef COSMO

#endif
#ifdef CESM1
   USE messy_main_data_bi,         ONLY: eps, q=>qm1, qm1
#endif
! ### add new BMs here

    IMPLICIT NONE
    INTRINSIC :: TRIM, ADJUSTL

    ! I/O
    CHARACTER(LEN=*),           INTENT(IN) :: vname ! input variable name
    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: dat   ! data

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'init_pvar'
    REAL(DP),         PARAMETER :: scvmr = M_air/M_H2O
    LOGICAL, SAVE               :: lfirst = .TRUE.
    INTEGER                     :: i
    CHARACTER(LEN=4)            :: nstr = ''
    CHARACTER(LEN=4)            :: sstr = ''
    LOGICAL                     :: l_cycle

    ! NO RE_INITIALISTION AFTER RESTART (PROGNOSTIC VARIABLES NEED TO BE
    ! IN RESTART-FILES)
    IF (lresume) RETURN

    IF (TRIM(bmstr) == '') THEN
       IF (lfirst) THEN
          CALL warning_bi(&
               'prognostic variable initialisation not'//&
               &' implemented for this basemodel' &
               ,substr)
       END IF
       lfirst = .FALSE.
       RETURN
    ENDIF

    DO i=1, N_INI_MAX

       CALL int2str(nstr, i, ' ', '*')

       ! CHECK BASEMODEL INFO
       l_cycle = .FALSE.
       SELECT CASE(TRIM(ADJUSTL(ini_pvar(i)%basemodel)))
       CASE('ECHAM5','COSMO','CESM1') ! ### add new BMs here
          ! OK
       CASE('')
          IF (lfirst) THEN
             CALL warning_bi(&
                  'empty basemodel name at #'//TRIM(nstr), substr)
          END IF
          l_cycle = .TRUE.
       CASE DEFAULT
          CALL error_bi(&
               'unknown basemodel name at #'//TRIM(nstr), substr)
       END SELECT
       IF (l_cycle) CYCLE

       l_cycle = .FALSE.
       ! CHECK PROGNOSTIC VARIABLE NAME
       IF (TRIM(ADJUSTL(ini_pvar(i)%pvar)) == '') THEN
          IF (lfirst) THEN
             CALL warning_bi(&
                  'empty prognostic variable name at #'//TRIM(nstr), substr)
          END IF
          l_cycle = .TRUE.
       ENDIF
       IF (l_cycle) CYCLE

       ! CHECK INPUT VARIABLE NAME
       IF (TRIM(ADJUSTL(ini_pvar(i)%inp)) == '') THEN
          CALL error_bi(&
               'empty prognostic variable name at #'//TRIM(nstr), substr)
       ENDIF
       !
       ! NOTHING ELSE TODO, IF INPUT VARIABLE NAME DOES NOT MATCH
       IF (TRIM(ADJUSTL(ini_pvar(i)%inp)) /= TRIM(ADJUSTL(vname))) CYCLE

       l_cycle = .FALSE.
       ! CHECK BASEMODEL
       IF (TRIM(ADJUSTL(ini_pvar(i)%basemodel)) /= TRIM(bmstr)) THEN
          CALL warning_bi(&
               'initialsiation of prognostic variable '//&
               &TRIM(ini_pvar(i)%pvar)//' at #'//TRIM(nstr)//&
               &'matches input variable, but not for basemodel '//TRIM(bmstr) &
               , substr)
          l_cycle = .TRUE.
       END IF
       IF (l_cycle) CYCLE

       CALL int2str(sstr, ini_pvar(i)%switch, ' ', '*')
       CALL info_bi('initialising prognostic variable '//&
            &TRIM(ini_pvar(i)%pvar)//' with'//&
            &' input variable '//TRIM(ini_pvar(i)%inp)//&
            &' (# '//TRIM(nstr)//', method '//TRIM(sstr)//')' ,substr)
       SELECT CASE(TRIM(ADJUSTL(ini_pvar(i)%pvar)))

#if defined(ECHAM5) || defined(CESM1)
          CASE('q')
             SELECT CASE(ini_pvar(i)%switch)
                CASE(0) ! input is already in [kg/kg]
                   q(:,:,:)   = dat(:,:,:)
                   qm1(:,:,:) = (1._DP - eps) * q(:,:,:)
                CASE(1) ! convert from [mol/mol] to [kg/kg]
                   q(:,:,:)   = dat(:,:,:)/(scvmr+dat(:,:,:))
                   qm1(:,:,:) = (1._DP - eps) * q(:,:,:)
                CASE DEFAULT
                   CALL error_bi('unknown method',substr)
             END SELECT
#if defined(CESM1)
             ! tell CESM1 basemodel that qm1 is already initialised
             ! USE from messy_main_tracer_mem_bi.f90
             qm1_init = .TRUE.
#endif
#endif

#ifdef COSMO

#endif

! ### add new BMS here

          CASE DEFAULT
             CALL error_bi(&
                  'unknown prognostic variable '//TRIM(ini_pvar(i)%pvar)//&
                  ' at #'//TRIM(nstr)//&
                  &' (or not yet implemented for basemodel'//&
                  &TRIM(bmstr)//')'&
                  , substr)

       END SELECT

    END DO

    lfirst = .FALSE.

  END SUBROUTINE init_pvar
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_read_nml_cpl(status, iou)

    ! TRACER MODULE ROUTINE (CORE)
    !
    ! READ TRACER NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Jul 2003

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

#if defined(ECHAM5)
!!#D attila +
    USE messy_attila_tools_e5, ONLY: LG2GP_SUM, LG2GP_AVE, LG2GP_STD &
                                   , LG2GP_AVEGT0
!!#D attila -
!!#D clams +
    USE messy_clams_tools_e5, ONLY: CL2GP_SUM
!!#D clams -
#endif

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CPL/ l_tracer_init, l_tracer_initfromrestart &
         , l_conv_lg2gp &
         , i_conv_lg2gp_mode, r_conv_lg2gp_fill, l_conv_lg2gp_mc &
         , ini_pvar     &
         , n_trcr_block &
         , l_conv_cl2gp, i_conv_cl2gp_mode, r_conv_cl2gp_fill

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_read_nml_cpl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

!!#D attila +
    IF (l_conv_lg2gp) THEN
       WRITE(*,*) '  CONVERSION OF LG TO GP       : ON '
       WRITE(*,*) '                    MODE       : ', i_conv_lg2gp_mode
#if defined(ECHAM5)
       SELECT CASE (i_conv_lg2gp_mode)
          CASE(LG2GP_SUM)
             WRITE(*,*) '                    ===>       : LG2GP_SUM'
          CASE(LG2GP_AVE)
             WRITE(*,*) '                    ===>       : LG2GP_AVE'
          CASE(LG2GP_STD)
             WRITE(*,*) '                    ===>       : LG2GP_STD'
          CASE(LG2GP_AVEGT0)
             WRITE(*,*) '                    ===>       : LG2GP_AVEGT0'
          CASE DEFAULT
             WRITE(*,*) '                    ===>       : UNKNOWN (ERROR)'
             RETURN
       END SELECT
#endif
       WRITE(*,*) '                    MODE       : ',i_conv_lg2gp_mode
       WRITE(*,*) '              FILL VALUE       : ',r_conv_lg2gp_fill
       IF (l_conv_lg2gp_mc) THEN
          WRITE(*,*) '       MASS CONSERVATION       : ON '
       ELSE
          WRITE(*,*) '       MASS CONSERVATION       : OFF '
       END IF
    ELSE
       WRITE(*,*) '  CONVERSION OF LG TO GP       : OFF'
    END IF
!!#D attila -

!!#D clams +
    IF (l_conv_cl2gp) THEN
       WRITE(*,*) '  CONVERSION OF CL TO GP       : ON '
       WRITE(*,*) '                    MODE       : ', i_conv_cl2gp_mode
#if defined(ECHAM5)
       SELECT CASE (i_conv_cl2gp_mode)
          CASE(CL2GP_SUM)
             WRITE(*,*) '                    ===>       : CL2GP_SUM'
          CASE(LG2GP_AVE)
             WRITE(*,*) '                    ===>       : CL2GP_AVE'
          CASE(LG2GP_STD)
             WRITE(*,*) '                    ===>       : CL2GP_STD'
          CASE(LG2GP_AVEGT0)
             WRITE(*,*) '                    ===>       : CL2GP_AVEGT0'
          CASE DEFAULT
             WRITE(*,*) '                    ===>       : UNKNOWN (ERROR)'
             RETURN
       END SELECT
#endif
       WRITE(*,*) '              FILL VALUE       : ',r_conv_cl2gp_fill
     ELSE
       WRITE(*,*) '  CONVERSION OF CL TO GP       : OFF'
     ENDIF
!!#D clams -

#ifdef COSMO
    WRITE(*,*) '   NUMBER OF TRACER BLOCKS:  ', n_trcr_block
#endif

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE main_tracer_read_nml_cpl
  ! -------------------------------------------------------------------

! **********************************************************************+
END MODULE messy_main_tracer_bi
! **********************************************************************+
