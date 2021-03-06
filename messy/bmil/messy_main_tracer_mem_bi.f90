! **********************************************************************+
MODULE messy_main_tracer_mem_bi
! **********************************************************************+

  ! MESSy
  USE messy_main_tracer, ONLY: dp, t_trinfo, t_trinfo_tp, t_trinfo_list &
       , modstr, ON, OFF, I_FORCE_COL, I_ADVECT, I_VDIFF, I_INTEGRATE, I_MIX &
       , I_HDIFF, I_RELAX ! um_ak_20081030 for COSMO
!#ifdef ICON
  USE messy_main_tracer, ONLY: STRLEN_TRSET
!#endif

  ! SPECIAL FOR LAGRANGIAN TRACERS / CHANNEL OBJECTS
!!#D attila +
#ifdef ECHAM5
  USE messy_attila_mem,  ONLY: NGCELL, NCELL
#endif
!!#D attila -

!!#D clams +
#ifdef ECHAM5
 USE messy_clams_global,  ONLY: dnparts_max, dnparts
#endif
!!#D clams -

  USE messy_main_tools,  ONLY: PTR_4D_ARRAY  ! um_ak_20081030

  IMPLICIT NONE
  !PUBLIC is already default
  PRIVATE :: dp
  ! op_bk_20130820+
#if defined(ICON)
  PUBLIC
#endif
  ! op_bk_20130820-
  SAVE

! op_pj_20120120+ moved to this position from messy_main_tracer_bi.f90
  ! GLOBAL SWITCHES FOR TRACER SETS
  LOGICAL :: L_GP   = .FALSE. ! GLOBAL SWITCH   ! mz_ab_20100226
  LOGICAL :: L_LG   = .FALSE. ! GLOBAL SWITCH
  LOGICAL :: L_OM   = .FALSE. ! GLOBAL SWITCH   ! mz_ap_20071023
! op_pj_20120120-
  LOGICAL :: L_CL   = .FALSE. ! GLOBAL SWITCH   ! ju_ec_20180611

#ifdef CESM1
  ! mz_ab_20171005+
  LOGICAL :: qm1_init = .FALSE.
  ! mz_ab_20171005-
#endif

  ! TRACER SET: GRIDPOINT ##################################################
  !
#if ! defined(ICON)
  ! NAME OF TRACER SET
  ! ub_ak_20190508+
  !CHARACTER(LEN=*), PARAMETER             :: GPTRSTR = 'gp'
  CHARACTER(LEN=STRLEN_TRSET), PARAMETER   :: GPTRSTR = 'gp'
  ! ub_ak_20190508-
  CHARACTER(LEN=*), PARAMETER              :: gp_channel = modstr//'_gp'
  !
#else
  !
  ! NAME OF TRACER SET and corresponding CHANNEL (for each patch)
  ! GPTRSTR    => L_GPTRSTR(dom_id)
  ! gp_channel => l_gp_channel(dom_id)
  CHARACTER(LEN=STRLEN_TRSET), POINTER                 :: GPTRSTR
!!$  CHARACTER(LEN=STRLEN_TRSET), POINTER                 :: gp_channel
  CHARACTER(LEN=*), PARAMETER              :: gp_channel = modstr//'_gp'
  CHARACTER(LEN=STRLEN_TRSET), DIMENSION(:), POINTER   :: L_GPTRSTR    => NULL()
!!$  CHARACTER(LEN=STRLEN_TRSET), DIMENSION(:), POINTER   :: l_gp_channel => NULL()
  INTEGER,                     DIMENSION(:), POINTER   :: l_ntrac_icon => NULL()
  INTEGER,                                   POINTER   :: ntrac_icon => NULL()
  !
#endif
  !
  ! POINTER TO TRACER META INFORMATION
  TYPE(t_trinfo_list),             POINTER :: trlist_gp => NULL()
  !
  ! ARRAY OF TRACER INFO STRUCTS
  TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti_gp => NULL()
  !
  ! NUMBER OF TRACERS
  INTEGER                                  :: ntrac_gp = 0
  !
  ! POINTER TO TRACER FIELDS
  ! - FOR INTERNAL MEMORY MANAGEMENT
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxt   => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxtte => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxtm1 => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxtf  => NULL()
  !
  ! - FOR OUTSIDE 'LOCAL LOOP'
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xt    => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtte  => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtm1  => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtf   => NULL()
  !
  ! - FOR INSIDE 'LOCAL LOOP'
  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxt   => NULL()
  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxtte => NULL()
  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxtm1 => NULL()
  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxtf  => NULL()
  !
  ! #######################################################################

  ! ju_ec_20180621+
  ! TRACER SET: LAGRANGE (CLaMS) ##########################################
  !
  ! NAME OF TRACER SET
  CHARACTER(LEN=*), PARAMETER             :: CLTRSTR = 'cl'
  CHARACTER(LEN=*), PARAMETER             :: cl_channel = modstr//'_cl'
  !
  ! POINTER TO TRACER META INFORAMTION
  TYPE(t_trinfo_list),             POINTER :: trlist_cl => NULL()
  !
  ! ARRAY OF TRACER INFO STRUCTS
  TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti_cl => NULL()
  !
  ! NUMBER OF TRACERS
  INTEGER                                  :: ntrac_cl = 0
  !
  ! POINTER TO TRACER FIELDS
  ! - FOR INTERNAL MEMORY MANAGEMENT
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER   :: pxt_c   => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER   :: pxtte_c => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER   :: pxtm1_c => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER   :: pxmem_c => NULL()
  ! - FOR OUTSIDE 'LOCAL LOOP'
  REAL(DP), DIMENSION(:,:,:,:),   POINTER   :: xt_c    => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER   :: xtte_c  => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER   :: xtm1_c  => NULL()
  ! - FOR INSIDE 'LOCAL LOOP'
  REAL(DP), DIMENSION(:,:),       POINTER   :: qxt_c   => NULL()
  REAL(DP), DIMENSION(:,:),       POINTER   :: qxtte_c => NULL()
  REAL(DP), DIMENSION(:,:),       POINTER   :: qxtm1_c => NULL()

  !
  ! FOR CL->GP CONVERSION
  CHARACTER(LEN=*), PARAMETER              :: CLGPTRSTR    = 'clgp'
  CHARACTER(LEN=*), PARAMETER              :: clgp_channel = modstr//'_clgp'
  TYPE(t_trinfo_list),             POINTER :: trlist_clgp => NULL()
  TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti_clgp     => NULL()
  INTEGER                                  :: ntrac_clgp  = 0
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER  :: pxt_clgp    => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER  :: xt_clgp     => NULL()
  REAL(DP), DIMENSION(:,:,:),     POINTER  :: qxt_clgp    => NULL()
  !
  ! #######################################################################
  ! ju_ec_20180621-


  ! TRACER SET: LAGRANGE (ATTILA) #########################################
!!#D attila +
#ifndef ECHAM5
!!#D attila -
  INTEGER, PARAMETER :: NCELL  = 0
  INTEGER, PARAMETER :: NGCELL = 0
!!#D attila +
#endif
!!#D attila -
  !
  ! NAME OF TRACER SET
  ! op_k_20190709+
  ! bugfix (g95 error: LGTRSTR and GPTRSTR need to have same length (see DDEP) )
  !CHARACTER(LEN=*), PARAMETER              :: LGTRSTR = 'lg'
  CHARACTER(LEN=STRLEN_TRSET), PARAMETER   :: LGTRSTR = 'lg'
  ! op_k_20190709-
  CHARACTER(LEN=*), PARAMETER              :: lg_channel = modstr//'_lg'
  !
  ! POINTER TO TRACER META INFORAMTION
  TYPE(t_trinfo_list),             POINTER :: trlist_lg => NULL()
  !
  ! ARRAY OF TRACER INFO STRUCTS
  TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti_lg => NULL()
  !
  ! NUMBER OF TRACERS
  INTEGER                                  :: ntrac_lg = 0
  !
  ! POINTER TO TRACER FIELDS
  ! - FOR INTERNAL MEMORY MANAGEMENT
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER   :: pxt_a   => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER   :: pxtte_a => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER   :: pxtm1_a => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER   :: pxtf_a  => NULL()
  ! - FOR OUTSIDE 'LOCAL LOOP'
  REAL(DP), DIMENSION(:,:,:,:),   POINTER   :: xt_a    => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER   :: xtte_a  => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER   :: xtm1_a  => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER   :: xtf_a   => NULL()
  ! - FOR INSIDE 'LOCAL LOOP'
  REAL(DP), DIMENSION(:,:),       POINTER   :: qxt_a   => NULL()
  REAL(DP), DIMENSION(:,:),       POINTER   :: qxtte_a => NULL()
  REAL(DP), DIMENSION(:,:),       POINTER   :: qxtm1_a => NULL()
  REAL(DP), DIMENSION(:,:),       POINTER   :: qxtf_a  => NULL()
  !
  ! FOR LG->GP CONVERSION
  CHARACTER(LEN=*), PARAMETER              :: LGGPTRSTR    = 'lggp'
  CHARACTER(LEN=*), PARAMETER              :: lggp_channel = modstr//'_lggp'
  TYPE(t_trinfo_list),             POINTER :: trlist_lggp => NULL()
  TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti_lggp     => NULL()
  INTEGER                                  :: ntrac_lggp  = 0
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER  :: pxt_lggp    => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER  :: xt_lggp     => NULL()
  REAL(DP), DIMENSION(:,:,:),     POINTER  :: qxt_lggp    => NULL()
  !
  ! DECOMPOSITION
  ! -> NUMBER OF CELLS ON THIS PE
  INTEGER                                   :: NLCELL    = 0
  ! -> CELLS ON THIS PE ?
  LOGICAL                                   :: LG_ACTIVE = .false.
  !
  INTEGER :: number_mix = 0
  ! #######################################################################

! mz_ap_20071023+
  ! TRACER SET: OCEANIC (MPIOM) #########################################
  !
  ! NAME OF TRACER SET
  CHARACTER(LEN=*), PARAMETER              :: OMTRSTR = 'om'
  CHARACTER(LEN=*), PARAMETER              :: om_channel = modstr//'_om'
  !
  ! POINTER TO TRACER META INFORAMTION
  TYPE(t_trinfo_list),             POINTER :: trlist_om => NULL()
  !
  ! ARRAY OF TRACER INFO STRUCTS
  TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti_om => NULL()
  !
  ! NUMBER OF TRACERS
  INTEGER                                  :: ntrac_om = 0
  !
  ! POINTER TO TRACER FIELDS
  ! - FOR INTERNAL MEMORY MANAGEMENT
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER   :: pxt_om   => NULL()
  ! - FOR USE IN SUBMODELS
  REAL(DP), DIMENSION(:,:,:,:),   POINTER   :: xt_om    => NULL()
  !
  ! #######################################################################
! mz_ap_20071023-

! um_ak_20080731+
  ! DEFINE BOUNDARY DATA VARIABLES FOR REGIONAL MODELS
  REAL(DP), DIMENSION(:,:,:,:,:),    POINTER :: xt_bd    => NULL()
  ! um_ak_20130521+
  ! FOR USE in VARTAB in COSMO/MESSy
  ! ACCESS of TIME LEVELS NECESSARY
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: xt_tl => NULL()
  ! um_ak_20130521-
  ! ub_ak_20170705+
  ! TRACER BLOCK STRUCTURE (required for COSMO only)
  ! NUMBER OF Blocked TRACERS
  INTEGER                                  :: ntracblck_gp = 0
  REAL(DP), DIMENSION(:,:,:,:),   POINTER  :: xtblck       => NULL()
  REAL(DP), DIMENSION(:,:,:),     POINTER  :: xtteblck     => NULL()
  ! ub_ak_20170705-
   ! #######################################################################
! um_ak_20080731-

! **********************************************************************+
END MODULE messy_main_tracer_mem_bi
! **********************************************************************+
