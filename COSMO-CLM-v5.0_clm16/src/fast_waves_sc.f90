!+ Module containing procedures to compute fast wave terms for RK-scheme
!------------------------------------------------------------------------------

MODULE fast_waves_sc

!------------------------------------------------------------------------------
!
! Description:
!   This module contains procedure(s) (mainly for the use in src_runge_kutta)
!   to update the prognostic variables
!   u, v, w, pp and T by the fast processes sound and gravity wave expansion
!   artificial divergence damping, and addition of the (constant) 
!   slow tendencies.
!
!   New aspects of this fast waves solver are
!   - correct weightings for all vertical derivatives and vertical averages
!   - strong conservation form of the divergence operator
!   - optional: 3D isotropic divergence filtering  
!   - optional: Mahrer (1984) discretization of horizontal pressure gradients
!
! Current Code Owner: DWD, Michael Baldauf
!  phone:  +49  69  8062 2733
!  fax:    +49  69  8062 3721
!  email:  michael.baldauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_24        2012/06/22 Michael Baldauf
!  Initial release
! V4_25        2012/09/28 Michael Baldauf, Hans-Juergen Panitz
!  Changed allowed values for itype_bbc_w: new nomenclature now requires itype_bbc_w > 100 (MB)
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep (HJP)
! V4_27        2013/03/19 Michael Baldauf
!  Use new namelist switch divdamp_slope
!  Consolidations with fast_waves_rk:
!   - optional use of the upper boundary damping layer by Klemp et al. (2008)
!   - lw_freeslip can be used
!   - new version of the radiative lateral boundary condition
! V4_28        2013/07/12 Michael Baldauf, Ulrich Schaettler
!  Improved 3D divergence damping by correct calling of SR 'l_calc_lhs_at_1st_RKstep'
!   and setting appropriate rhs (dependent on switches l_no_div_damping_in_1st_step, ...).
!  Adaptation of the dynamical bottom boundary condition to the newer formulation
!   of the buoyancy term and small bug fix in the use of dw_dt_slow.
!   (changes in results only, if ldyn_bbc=.TRUE. or 3D div.-damping is used.
!    which must be activated by an internal switch l_3D_div_damping=.TRUE.)
!  Use vcoord%kflat from new module vgrid_refatm_utils (US)
! V4_29        2013/10/04 Ulrich Schaettler
!  Replaced vcoord_in by use of vcoord
!  Initialize izerror in all subroutines
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!============================================================================
!
! Declarations:
!
! Modules used:
!

USE data_parameters , ONLY :   &
  ireals,    & ! KIND-type parameters for real variables
  iintegers    ! KIND-type parameter for standard integer variables

!----------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &
  
  ! 2. horizontal and vertical sizes of the fields and related variables
  ! --------------------------------------------------------------------
  ie,           & ! number of grid points in zonal direction
  je,           & ! number of grid points in meridional direction
  ke,           & ! number of grid points in vertical direction
  ke1,          & ! ke+1
  ke_rd,        & ! lowest level with Rayleigh-damping
  
  ! 3. start- and end-indices for the computations in the horizontal layers
  ! -----------------------------------------------------------------------
  !    These variables give the start- and the end-indices of the
  !    forecast for the prognostic variables in a horizontal layer.
  !    Note, that the indices for the wind-speeds u and v differ from
  !    the other ones because of the use of the staggered Arakawa-B-grid.
  !
  !    zonal direction
  istart,       & ! start index for the forecast of w, t, qv, qc and pp
  iend,         & ! end index for the forecast of w, t, qv, qc and pp
  istartu,      & ! start index for the forecast of u
  iendu,        & ! end index for the forecast of u
  istartv,      & ! start index for the forecast of v
  iendv,        & ! end index for the forecast of v
  istartpar,    & ! start index for computations in the parallel program
  iendpar,      & ! end index for computations in the parallel program
  
  !    meridional direction
  jstart,       & ! start index for the forecast of w, t, qv, qc and pp
  jend,         & ! end index for the forecast of w, t, qv, qc and pp
  jstartu,      & ! start index for the forecast of u
  jendu,        & ! start index for the forecast of u
  jstartv,      & ! start index for the forecast of v
  jendv,        & ! end index for the forecast of v
  jstartpar,    & ! start index for computations in the parallel program
  jendpar,      & ! end index for computations in the parallel program
  
  ! 4. constants for the horizontal rotated grid and related variables
  ! ------------------------------------------------------------------
  dlon,         & ! in degrees
  dlat,         & ! in degrees
  eddlon,       & ! 1 / dlon (in 1/Bogenmass)
  eddlat,       & ! 1 / dlat (in 1/Bogenmass)

  ! 5. variables for the time discretization and related variables
  ! --------------------------------------------------------------
  
  betasw,       & ! beta-variable for treatment of soundwaves    in w
  betagw,       & ! beta-variable for treatment of gravity-waves in w
  beta2sw,      & ! beta-variable for treatment of soundwaves    in p*, T*
  beta2gw         ! beta-variable for treatment of gravity-waves in p*, T*

! end of data_modelconfig

!----------------------------------------------------------------------------

USE data_fields     , ONLY :   &
  rho0       ,    & ! reference density at the full model levels    (kg/m3)
  p0         ,    & ! reference pressure    at main levels          ( Pa  )
  t0         ,    & ! reference temperature at main levels          ( K   )
  dt0dz      ,    & ! temperature gradient of reference atmosphere  ( K/m )
  hhl        ,    & ! geometrical height of half model levels       ( m   )
  crlat      ,    & ! cosine of transformed latitude
  acrlat     ,    & ! 1 / ( crlat * radius of the earth )           ( 1/m )
  rdcoef            ! Rayleigh damping coefficients 

USE grid_metrics_utilities,  ONLY:  &
  sqrtg_r_s  ,    & ! reciprocal square root of G at skalar points  ( 1/m )
  sqrtg_r_u  ,    & ! reciprocal square root of G at u points       ( 1/m )
  sqrtg_r_v  ,    & ! reciprocal square root of G at v points       ( 1/m )
  sqrtg_r_w  ,    & ! reciprocal square root of G at w points       ( 1/m )
  dzeta_dlam ,    &
  dzeta_dphi ,    &
  wgtfac     ,    &
  wgtfac_u   ,    &
  wgtfac_v   ,    &
  wgt_one_sided


! end of data_fields

!----------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
  num_compute,     & ! number of compute PEs
  nboundlines,     & ! number of boundary lines of the domain for which
                     ! no forecast is computed = overlapping boundary
  ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
  ncomm_type,      & ! type of communication
  my_cart_id,      & ! rank of this subdomain in the cartesian communicator
  my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
  icomm_cart,      & ! communicator for the virtual cartesian topology
  imp_reals,       & ! determines the correct REAL type used in the model
                     ! for MPI
  imp_integers,    & ! determines the correct INTEGER type used in the
                     ! model for MPI
  nexch_tag,       & ! tag to be used for MPI boundary exchange
                     !  (in calls to exchg_boundaries)
  sendbuf,         & ! sending buffer for boundary exchange:
                     ! 1-4 are used for sending, 5-8 are used for receiving
  isendbuflen        ! length of one column of sendbuf

!----------------------------------------------------------------------------

USE parallel_utilities, ONLY :  &
  global_values

USE data_runcontrol, ONLY :   &
  
  ! 4. controlling the dynamics
  ! ---------------------------
  divdamp_slope, & ! exceed the theoretical slope stability criterion of
                   ! the divergence damping (only for itype_fast_waves=2)

  ! 7. additional control variables
  ! -------------------------------
  lradlbc,      & ! if lartif_data=.TRUE.: radiative lateral boundary cond.
                  !               .FALSE.:   or with Davies conditions
  l2dim,        & ! if lartif_data=.TRUE.: 2dimensional model version
  lperi_x,      & ! if lartif_data=.TRUE.: periodic boundary cond. in x-dir.
                  !               .FALSE.:   or with Davies conditions
  lperi_y,      & ! if lartif_data=.TRUE.: periodic boundary cond. in y-dir.
                  !               .FALSE.:   or with Davies conditions
  itype_bbc_w,  & ! Bottom boundary condition for vertical wind
                  ! =0/1: RK-like method following iadv_order
                  ! =2/3: differencing following iadv_order without RK stepping
                  ! =4/5: Fourth-order centered differences
                  ! 0/2/4: include extrapolation of horizontal wind to sfc
                  ! 1/3/5: no extrapolation of horizontal wind to sfc
                  ! alternative nomenclature: if >= 100, i.e. =1ed then
                  !   e = order of extrapolation (0,1,2)
                  !   d = order of dh/dx, dh/dy (2,3,4)
  lw_freeslip,  & ! if .TRUE.: with free slip lateral boundary condition and
                  ! if .FALSE. specified lateral boundary values for w
  ldyn_bbc,     & ! if .TRUE., dynamical bottom boundary condition
  lspubc ,      & ! if .TRUE., use Rayleigh damping in the upper levels
  itype_spubc,  & ! type of Rayleigh damping in the upper levels
  irunge_kutta, &
  idbg_level,   &
  lprintdeb_all,&
  ldebug_dyn

!-----------------------------------------------------------------------------

USE data_constants  , ONLY :   &
  
  ! 2. physical constants and related variables
  ! -------------------------------------------
  r_d,          & ! gas constant for dry air
  cp_d,         &
  rvd_m_o,      & ! r_v/r_d - 1; hier fuer Subr. calrho
  gamma,        & ! 1 / (1 - rdocp)   ( = cp_d/cv_d)
  r_earth,      & ! mean radius of the earth (m)
  g,            & ! acceleration due to gravity (m/s^2)
  pi              ! circle constant

!-----------------------------------------------------------------------------

USE environment     , ONLY :  &
  exchg_boundaries,    &
  model_abort

USE vgrid_refatm_utils, ONLY: &
  vcoord          ! definitions for vertical grid (here: kflat)

USE data_runcontrol, ONLY: ltime

!=============================================================================

IMPLICIT NONE

!=============================================================================

! private module variables:


! --- fields for metrics:

REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: dz_dlam (:,:,:)
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: dz_dphi (:,:,:)

! 4th order approximation (only at the bottom k=ke+1):
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: dz_dlam_4ord (:,:)
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: dz_dphi_4ord (:,:)

REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: dz_lam_u_half (:,:,:)
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: dz_phi_v_half (:,:,:)

REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: sqrtgtot_r       (:,:,:)


! --- switches, fields for the linear equation system:

REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: lgs_A   (:,:,:)
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: lgs_B   (:,:,:)
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: lgs_C   (:,:,:)
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: lgs_rhs (:,:,:)

LOGICAL :: l_calc_lhs_at_1st_RKstep
! if .TRUE. then call subroutine 'lhs_of_tridiag_system_for_w'
! at the beginning of a large RK-time step
! this switch depends on l_small_pert_in_pT !

LOGICAL, PRIVATE :: l_use_tridag
! =.TRUE. : solve the tridiagonal system via the numerical recipes 
!           subroutine 'tridag'
! =.FALSE.: use an explicitely programmed LU-decomposition


! --- Crank-Nicholson weighting parameters:

REAL (KIND=ireals), PRIVATE ::  &
     beta_s_1_impl, beta_s_2_impl, beta_s_3_impl, &
     beta_s_4_impl, beta_s_5_impl, beta_s_6_impl, &
     beta_s_7_impl, beta_s_8_impl, beta_s_9_impl

REAL (KIND=ireals), PRIVATE ::  &
     beta_b_1_impl, beta_b_2_impl, beta_b_3_impl, beta_b_4_impl

REAL (KIND=ireals), PRIVATE ::  &
     beta_d_1_impl, beta_d_2_impl, beta_d_3_impl,  &
     beta_d_4_impl, beta_d_5_impl, beta_d_6_impl

REAL (KIND=ireals), PRIVATE ::  &
     beta_s_1_expl, beta_s_2_expl, beta_s_3_expl,  &
     beta_s_4_expl, beta_s_5_expl, beta_s_6_expl,  &
     beta_s_7_expl, beta_s_8_expl, beta_s_9_expl

REAL (KIND=ireals), PRIVATE ::  &
     beta_b_1_expl, beta_b_2_expl, beta_b_3_expl, beta_b_4_expl

REAL (KIND=ireals), PRIVATE ::  &
     beta_d_1_expl, beta_d_2_expl, beta_d_3_expl,    &
     beta_d_4_expl, beta_d_5_expl, beta_d_6_expl

INTEGER (KIND=iintegers), PRIVATE :: set_beta_s_8_9_to


! --- thermodynamics constants, fields:

REAL (KIND=ireals), PRIVATE :: cv_d, R_by_cv

REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: T_tot            (:,:,:)
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: gamma_p_tot      (:,:,:)  ! = c_p / c_v * p_tot
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: rho              (:,:,:)
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: rho_inv          (:,:,:)
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: frac_1_rhodx     (:,:,:)
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: frac_1_rhody     (:,:,:)
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: rho_inv_k_avg    (:,:,:)
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: p0T_pT0_k_avg    (:,:,:)  ! = (p0 * T ) / (p * T0)
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: p0_pT0_k_avg     (:,:,:)  ! = p0 / (p * T0)
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: p_inv_k_avg      (:,:,:)  ! = 1/p

LOGICAL, PRIVATE :: l_small_pert_in_pT 
! = .TRUE.: assume small perturbations in the environment,
! i.e. the calculation of all coefficients
! which are functions of T or p are done outside of the time integration loop.
! consequence: the lhs of the tridiagonal system must be calculated only once!
! In very fast varying environments (e.g. tornado simulations??)
! it could be advisable to set it .FALSE. 


! --- Variables, fields for divergence damping:

LOGICAL          :: l_3D_div_damping
LOGICAL, PRIVATE :: l_no_div_damping_in_1st_step
LOGICAL, PRIVATE :: limit_alpha_div_by_slope_stability

REAL (KIND=ireals), ALLOCATABLE :: alpha_div_h      (:,:,:)
REAL (KIND=ireals) :: alpha_div_v_to_h


! --- Variables, fields for Mahrer (1984) discretization:

LOGICAL, PRIVATE :: l_hor_p_grad_Mahrer
LOGICAL, PRIVATE :: l_Mahrer_opti      !MB: sollte spaeter wieder rausfallen
   ! tolerance thickness for the weak extrapolation beneath the ground:
REAL   (KIND=ireals), PRIVATE :: delta_z_Mahrer
INTEGER(KIND=iintegers), PRIVATE :: k_index_same_level_Mahrer
LOGICAL, PRIVATE       :: is_on_same_level
LOGICAL, PRIVATE       :: is_search_finished

INTEGER(KIND=iintegers), ALLOCATABLE, PRIVATE :: p_hor_idx(:,:,:,:)
REAL   (KIND=ireals),    ALLOCATABLE, PRIVATE :: p_hor_wght(:,:,:,:)
INTEGER(KIND=iintegers), ALLOCATABLE, PRIVATE :: k_star_u(:,:)
INTEGER(KIND=iintegers), ALLOCATABLE, PRIVATE :: k_star_v(:,:)


! --- variables driving boundary conditions/treatment:

! Weighting factors for Rayleigh damping by Klemp et al. (2008) 
REAL (KIND=ireals), ALLOCATABLE, PRIVATE :: wrdfac(:)

LOGICAL, PRIVATE :: l_use_lateral_bound_values

INTEGER, PRIVATE :: i_bbc_w_extrapol_uv_order
INTEGER, PRIVATE :: i_bbc_w_dhdx_order

INTEGER (KIND=iintegers), PRIVATE :: itype_bound_for_Z_horiz

INTEGER, PRIVATE :: i_div_at_horiz_bound 
! =1:    extra treatment in k=1 and k=ke
! =2: no extra treatment in k=1 and k=ke

! --- other variables:

INTEGER, PRIVATE :: i_start_u
INTEGER, PRIVATE :: j_start_v
INTEGER, PRIVATE :: izdebug

!============================================================================

CONTAINS

!============================================================================

SUBROUTINE init_fast_waves_sc_1

  !--------------------------------------------------------------------------
  !
  ! Description:
  !   Initialization subroutine (1).
  !   It is called only once at the begin of the simulation.
  !
  !--------------------------------------------------------------------------

  INTEGER :: i, j, k
  LOGICAL :: is_div_damp_VI
  INTEGER :: istat

  REAL (KIND=ireals) :: zu, zv, zp1, zp2
  INTEGER            :: m

  REAL    (KIND=ireals),    ALLOCATABLE :: hml(:,:,:)  ! height of main levels
  INTEGER (KIND=iintegers), ALLOCATABLE :: k_star(:,:,:)

  CHARACTER(LEN=100) :: yzerrmsg
  INTEGER            :: izerror

  IF (ldebug_dyn) THEN
    IF (lprintdeb_all) THEN
      izdebug = idbg_level
    ELSE
      IF (my_cart_id == 0) THEN
        izdebug = idbg_level
      ELSE
        izdebug = 0
      ENDIF
    ENDIF
  ELSE
    izdebug = 0
  ENDIF

  IF (izdebug >= 10) THEN
    WRITE(*,*) "Subr. [init_fast_waves_sc_1]..."
  END IF

  ! some abbreviations:
  !-----------------------------------------------------------------

  cv_d    = cp_d - R_d
  R_by_cv = R_d / cv_d   ! = gamma - 1


  ! choose switches
  ! ----------------------------------------------------------------

  limit_alpha_div_by_slope_stability = .TRUE.
  l_3D_div_damping                   = .FALSE.
  l_no_div_damping_in_1st_step       = .TRUE.

  l_use_tridag                       = .TRUE.

  set_beta_s_8_9_to                  = 46

  itype_bound_for_Z_horiz            = 2

  l_hor_p_grad_Mahrer                = .FALSE.

  i_div_at_horiz_bound               = 2


  ! plausibility of lateral boundary treatment:
  
  ! periodic BC's have highest priority:
  IF ( lperi_x .OR. lperi_y ) THEN
    IF ( lradlbc .AND. (my_cart_id == 0) ) THEN
      WRITE(*,*) "WARNING: lradlbc is set to .FALSE.!"
    END IF
    lradlbc = .FALSE.

    IF (lw_freeslip .AND. (my_cart_id == 0) ) THEN
      WRITE(*,*) "WARNING: lw_freeslip is set to .FALSE.!"
    END IF
    lw_freeslip = .FALSE.
  END IF

  ! radlbc has higher priority than bound. values
  IF ( lradlbc ) THEN
    l_use_lateral_bound_values = .FALSE.
  ELSE
    l_use_lateral_bound_values = .TRUE.
  END IF

  is_div_damp_VI = .TRUE.

  l_small_pert_in_pT       = .TRUE.

  ! distinguish the following cases:
  ! 1.) 'fast' change in T and p can occur:
  !   --> l_small_pert_in_pT =.FALSE.
  !        --> l_calc_lhs_at_1st_RKstep=.FALSE.
  ! 2.) 'fast' change in T and p is not expected (default case):
  !   --> l_small_pert_in_pT =.TRUE.
  !    a.) if  i_runge_kutta=1:
  !        1.) if  l_3D_div_damping=.FALSE.
  !            --> l_calc_lhs_at_1st_RKstep=.TRUE.
  !        2.) if  l_3D_div_damping=.TRUE.
  !            a.) if  l_no_div_damping_in_1st_step=.TRUE.
  !                 --> l_calc_lhs_at_1st_RKstep=.FALSE.
  !            b.) if l_no_div_damping_in_1st_step=.FALSE.
  !        --> l_calc_lhs_at_1st_RKstep=.TRUE.
  !    b.) if  i_runge_kutta=2:
  !        --> l_calc_lhs_at_1st_RKstep=.FALSE.
  !

  IF ( .NOT. l_small_pert_in_pT ) THEN
    ! in this case, you can't calculate the coefficients
    ! only once per large timestep
    l_calc_lhs_at_1st_RKstep = .FALSE.
  ELSE
    IF ( irunge_kutta == 1 ) THEN
      IF ( .NOT. l_3D_div_damping ) THEN
        l_calc_lhs_at_1st_RKstep = .TRUE.  ! this is the normal case
      ELSE
        IF ( l_no_div_damping_in_1st_step ) THEN
          l_calc_lhs_at_1st_RKstep = .FALSE.
        ELSE
          l_calc_lhs_at_1st_RKstep = .TRUE.
        END IF
      END IF
    ELSE IF ( irunge_kutta == 2 ) THEN
      l_calc_lhs_at_1st_RKstep = .FALSE.
    ELSE
      l_calc_lhs_at_1st_RKstep = .FALSE.
    END IF
  END IF

  ! set weighting factors for the implicit (Crank-Nicholson) scheme:
  !----------------------------------------------------------------------

  ! betasw  is the Ikawa beta parameter for soundwaves    in w
  ! betagw            "                 for gravity-waves in w
  ! beta2sw           "                 for soundwaves    in p*, T*
  ! beta2gw           "                 for gravity-waves in p*, T*
  ! (=0: centered, =1 backward, =-1 forward)

  ! sound vertical implicit:
  beta_s_1_impl = 0.0_ireals   ! necessary for maximal tridiag. system for w
  beta_s_2_impl = ( 1.0_ireals + betasw  ) / 2.0_ireals
  beta_s_3_impl = 1.0_ireals   ! necessary for forward-backward scheme
  beta_s_4_impl = ( 1.0_ireals + beta2sw ) / 2.0_ireals
  beta_s_5_impl = 1.0_ireals   ! necessary for forward-backward scheme
  beta_s_6_impl = ( 1.0_ireals + beta2sw ) / 2.0_ireals
  beta_s_7_impl = 0.0_ireals   ! necessary for maximal tridiag. system for w

  !beta_s_8_impl = 1.0_ireals 
  !beta_s_9_impl = 1.0_ireals

  IF ( set_beta_s_8_9_to == 35 ) THEN
    IF ( my_cart_id == 0 ) THEN
      WRITE(*,*) "set  beta_s_8=beta_s_3, ... "
    ENDIF

    beta_s_8_impl = beta_s_3_impl
    beta_s_9_impl = beta_s_5_impl

  ELSE IF ( set_beta_s_8_9_to == 46 ) THEN
    IF ( my_cart_id == 0 ) THEN
      WRITE(*,*) "set  beta_s_8=beta_s_4, ... ", &
               "(recommended together with i_div_at_horiz_bound=1)"
    ENDIF
    beta_s_8_impl = beta_s_4_impl
    beta_s_9_impl = beta_s_6_impl
  ELSE
    yzerrmsg = "false value in set_beta_s_8_9_to"
    CALL model_abort (my_cart_id, 1234, yzerrmsg, "init_fast_waves_sc_1")
  END IF

  IF ( i_div_at_horiz_bound == 1 .AND. ( set_beta_s_8_9_to == 35 ) ) THEN
    IF ( my_cart_id == 0 ) THEN
      WRITE(*,*) "======================================================="
      WRITE(*,*) " WARNING: for i_div_at_horiz_bound==1 also "
      WRITE(*,*) "          set_beta_s_8_9_to == 46 must be set!"
      WRITE(*,*) "======================================================="
    ENDIF
  END IF

  ! buoyancy vertically implicit:
  beta_b_1_impl = ( 1.0_ireals + betagw  ) / 2.0_ireals
  beta_b_2_impl = ( 1.0_ireals + betagw  ) / 2.0_ireals
  beta_b_3_impl = ( 1.0_ireals + beta2gw ) / 2.0_ireals
  beta_b_4_impl = ( 1.0_ireals + beta2gw ) / 2.0_ireals

  IF ( is_div_damp_VI ) THEN
    IF ( my_cart_id == 0 ) THEN
      WRITE(*,*) "divergence damping vertically implicit"
    ENDIF
    beta_d_1_impl = 0.0_ireals   ! necessary for maximal tridiag. system for w
    beta_d_2_impl = 0.0_ireals   ! necessary for maximal tridiag. system for w
    beta_d_3_impl = 1.0_ireals 
    beta_d_4_impl = 1.0_ireals  
    beta_d_5_impl = 0.0_ireals   ! necessary for maximal tridiag. system for w
    !beta_d_6_impl = 1.0_ireals  
    beta_d_6_impl = beta_d_4_impl
  ELSE
    IF ( my_cart_id == 0 ) THEN
      WRITE(*,*) "divergence damping forward-backward"
    ENDIF
    beta_d_1_impl = 0.0_ireals   ! necessary for maximal tridiag. system for w
    beta_d_2_impl = 0.0_ireals   ! necessary for maximal tridiag. system for w
    beta_d_3_impl = 1.0_ireals 
    beta_d_4_impl = 0.0_ireals  
    beta_d_5_impl = 0.0_ireals   ! necessary for maximal tridiag. system for w
    !beta_d_6_impl = 1.0_ireals
    beta_d_6_impl = beta_d_4_impl
  END IF

  ! set the explicit counterparts:
  beta_s_1_expl = 1.0_ireals - beta_s_1_impl
  beta_s_2_expl = 1.0_ireals - beta_s_2_impl
  beta_s_3_expl = 1.0_ireals - beta_s_3_impl
  beta_s_4_expl = 1.0_ireals - beta_s_4_impl
  beta_s_5_expl = 1.0_ireals - beta_s_5_impl
  beta_s_6_expl = 1.0_ireals - beta_s_6_impl
  beta_s_7_expl = 1.0_ireals - beta_s_7_impl
  beta_s_8_expl = 1.0_ireals - beta_s_8_impl
  beta_s_9_expl = 1.0_ireals - beta_s_9_impl

  beta_b_1_expl = 1.0_ireals - beta_b_1_impl
  beta_b_2_expl = 1.0_ireals - beta_b_2_impl
  beta_b_3_expl = 1.0_ireals - beta_b_3_impl
  beta_b_4_expl = 1.0_ireals - beta_b_4_impl

  beta_d_1_expl = 1.0_ireals - beta_d_1_impl
  beta_d_2_expl = 1.0_ireals - beta_d_2_impl
  beta_d_3_expl = 1.0_ireals - beta_d_3_impl
  beta_d_4_expl = 1.0_ireals - beta_d_4_impl
  beta_d_5_expl = 1.0_ireals - beta_d_5_impl
  beta_d_6_expl = 1.0_ireals - beta_d_6_impl



  ! to be sure:
  vcoord%kflat = MAX( 1, vcoord%kflat )    ! upper half level (k=1) is flat in any case

  IF ( l_hor_p_grad_Mahrer ) THEN

    ! For each u- or v-gridpoint (i,j,k)
    ! estimate the vertical indices p_hor_idx(l,i,j,k)
    ! (l=1,2 for u, l=3,4 for v) 
    ! and linear interpolation weights p_hor_wght(l,i,j,k)
    ! of those 'scalar' gridpoints (=position of p') which lie
    ! around the point of the same height.
    ! Additionally k_star_u, k_star_v
    ! denote the lowest level, on which an extrapolation is advisable.
    ! (s. Zaengl (2012) MWR) 

    ALLOCATE( p_hor_idx (4,ie,je,ke), STAT=istat )
    ALLOCATE( p_hor_wght(4,ie,je,ke), STAT=istat )

    ALLOCATE( hml(ie,je,ke), STAT=istat )

    ALLOCATE( k_star(4,ie,je), STAT=istat )

    ALLOCATE( k_star_u(ie,je), STAT=istat )
    ALLOCATE( k_star_v(ie,je), STAT=istat )


    ! tolerance thickness for the weak extrapolation beneath the ground:
    delta_z_Mahrer = 20.0 ! meter (see Zaengl (2012) MWR)

    DO k=1, ke
      DO j=1, je
        DO i=1, ie
          hml(i,j,k) = 0.5_ireals * ( hhl(i,j,k) + hhl(i,j,k+1) )
        END DO
      END DO
    END DO

    k_star(:,:,:) = 1

    DO k=1, ke
      DO j=2, je-1
        DO i=2, ie-1

          zu = 0.5_ireals * ( hml(i,j,k) + hml(i+1,j,k) )

          ! p-position right (eastward) of u(i,j,k):

          IF ( zu >= hml(i+1,j,1) ) THEN
            ! extrapolate p' upwards:
            zp1 = hml(i+1,j,1)
            zp2 = hml(i+1,j,2)
            p_hor_idx(1,i,j,k) = 2   ! Spezialfall
            k_star(1,i,j) = MAX( k_star(1,i,j), k)

          ELSE IF ( zu >= hml(i+1,j,ke) ) THEN
            ! interpolation is possible:
            DO m=2, ke 
              IF ( ( hml(i+1,j,m)  <= zu ) .AND. (zu < hml(i+1,j,m-1)) ) THEN 
                zp1 = hml(i+1,j,m-1)
                zp2 = hml(i+1,j,m  )
                p_hor_idx(1,i,j,k) = m
                k_star(1,i,j) = MAX( k_star(1,i,j), k)
                EXIT
              END IF
            END DO

          ELSE IF ( zu >= hml(i+1,j,ke) - delta_z_Mahrer ) THEN
            ! slight extrapolation of p' downwards:
            zp1 = hml(i+1,j,ke-1)
            zp2 = hml(i+1,j,ke  )
            p_hor_idx(1,i,j,k) = ke+1  ! special case
            k_star(1,i,j) = MAX( k_star(1,i,j), k)

          ELSE
            ! extrapolation downwards is possible, but not advisable:
            zp1 = hml(i+1,j,ke-1)
            zp2 = hml(i+1,j,ke  )
            p_hor_idx(1,i,j,k) = ke+1  ! special case

          END IF
          p_hor_wght(1,i,j,k) = ( zu - zp2 ) / ( zp1 - zp2 ) 


          ! p-position left (westward) of u(i,j,k):

          IF ( zu >= hml(i,j,1) ) THEN
            ! extrapolate p' upwards:
            zp1 = hml(i,j,1)
            zp2 = hml(i,j,2)
            p_hor_idx(2,i,j,k) = 2   ! special case
            k_star(2,i,j) = MAX( k_star(2,i,j), k)

          ELSE IF ( zu >= hml(i,j,ke) ) THEN
            ! interpolation is possible:
            DO m=2, ke 
              IF ( ( hml(i,j,m)  <= zu ) .AND.  ( zu < hml(i,j,m-1) ) ) THEN 
                zp1 = hml(i,j,m-1)
                zp2 = hml(i,j,m  )
                p_hor_idx(2,i,j,k) = m
                k_star(2,i,j) = MAX( k_star(2,i,j), k)
                EXIT
              END IF
            END DO

          ELSE IF ( zu >= hml(i,j,ke) - delta_z_Mahrer ) THEN
            ! slight extrapolation of p' downwards:
            zp1 = hml(i,j,ke-1)
            zp2 = hml(i,j,ke  )
            p_hor_idx(2,i,j,k) = ke+1  ! Spezialfall
            k_star(2,i,j) = MAX( k_star(2,i,j), k)

          ELSE
            ! extrapolation downwards is possible, but not advisable:
            zp1 = hml(i,j,ke-1)
            zp2 = hml(i,j,ke  )
            p_hor_idx(2,i,j,k) = ke+1  ! special case

          END IF
          p_hor_wght(2,i,j,k) = ( zu - zp2 ) / ( zp1 - zp2 ) 


          zv = 0.5_ireals * ( hml(i,j,k) + hml(i,j+1,k) )

          ! p-position northward of v(i,j,k):
          IF ( zv >= hml(i,j+1,1) ) THEN
            ! extrapoliere p' nach oben:
            zp1 = hml(i,j+1,1)
            zp2 = hml(i,j+1,2)
            p_hor_idx(3,i,j,k) = 2   ! special case
            k_star(3,i,j) = MAX( k_star(3,i,j), k)

          ELSE IF ( zv >= hml(i,j+1,ke) ) THEN
            ! interpolation is possible:
            DO m=2, ke 
              IF ( ( hml(i,j+1,m) <= zv ) .AND.  (zv < hml(i,j+1,m-1)) ) THEN 
                zp1 = hml(i,j+1,m-1)
                zp2 = hml(i,j+1,m  )
                p_hor_idx(3,i,j,k) = m
                k_star(3,i,j) = MAX( k_star(3,i,j), k)
                EXIT
              END IF
            END DO

          ELSE IF ( zv >= hml(i,j+1,ke) - delta_z_Mahrer ) THEN
            ! slight extrapolation of p' downwards:
            zp1 = hml(i,j+1,ke-1)
            zp2 = hml(i,j+1,ke  )
            p_hor_idx(3,i,j,k) = ke+1  ! special case
            k_star(3,i,j) = MAX( k_star(3,i,j), k)

          ELSE
            ! extrapolation downwards is possible, but not advisable:
            zp1 = hml(i,j+1,ke-1)
            zp2 = hml(i,j+1,ke  )
            p_hor_idx(3,i,j,k) = ke+1  ! Spezialfall

          END IF
          p_hor_wght(3,i,j,k) = ( zv - zp2 ) / ( zp1 - zp2 ) 


          ! p-position southward of v(i,j,k):

          IF ( zv >= hml(i,j,1) ) THEN
            ! extrapoliere p' nach oben:
            zp1 = hml(i,j,1)
            zp2 = hml(i,j,2)
            p_hor_idx(4,i,j,k) = 2   ! special case
            k_star(4,i,j) = MAX( k_star(4,i,j), k)

          ELSE IF ( zv >= hml(i,j,ke) ) THEN
            ! interpolation is possible:
            DO m=2, ke 
              IF ( ( hml(i,j,m)  <= zv ) .AND.  (zv < hml(i,j,m-1) ) ) THEN 
                zp1 = hml(i,j,m-1)
                zp2 = hml(i,j,m  )
                p_hor_idx(4,i,j,k) = m
                k_star(4,i,j) = MAX( k_star(4,i,j), k)
                EXIT
              END IF
            END DO

          ELSE IF ( zv >= hml(i,j,ke) - delta_z_Mahrer ) THEN
            ! slight extrapolation of p' downwards:
            zp1 = hml(i,j,ke-1)
            zp2 = hml(i,j,ke  )
            p_hor_idx(4,i,j,k) = ke+1  ! special case
            k_star(4,i,j) = MAX( k_star(4,i,j), k)

          ELSE
            ! extrapolation downwards is possible, but not advisable:
            zp1 = hml(i,j,ke-1)
            zp2 = hml(i,j,ke  )
            p_hor_idx(4,i,j,k) = ke+1  ! special case

          END IF
          p_hor_wght(4,i,j,k) = ( zv - zp2 ) / ( zp1 - zp2 ) 


        END DO
      END DO
    END DO

    ! k_star = max. k, for which a horizontal gradient can be determined by 
    ! interpolation or (weak) extrapolation on both sides of a 
    ! u- or a v-grid point, respectively (Zaengl (2012) MWR)
    DO j=2, je-1
      DO i=2, ie-1
        k_star_u(i,j) = MIN( k_star(1,i,j), k_star(2,i,j) )
        k_star_v(i,j) = MIN( k_star(3,i,j), k_star(4,i,j) )
      END DO
    END DO

    ! Until which level k do all the pressure points are lying in the same
    ! layer as their appropriate u or v points (for the 
    ! horizontal pressure calculation)?

    k_index_same_level_Mahrer = vcoord%kflat-1
    is_search_finished = .FALSE.

    DO k=vcoord%kflat, ke

      DO j=2, je-1
        DO i=2, ie-1

          !is_on_same_level =                                                &
          !        ((p_hor_idx(1,i,j,k)==k) .OR. (p_hor_idx(1,i,j,k)==k+1) ) &
          !  .AND. ((p_hor_idx(2,i,j,k)==k) .OR. (p_hor_idx(2,i,j,k)==k+1) ) &
          !  .AND. ((p_hor_idx(3,i,j,k)==k) .OR. (p_hor_idx(3,i,j,k)==k+1) ) &
          !  .AND. ((p_hor_idx(4,i,j,k)==k) .OR. (p_hor_idx(4,i,j,k)==k+1) )

          zu = 0.5_ireals * ( hml(i,j,k) + hml(i+1,j,k) )
          zv = 0.5_ireals * ( hml(i,j,k) + hml(i,j+1,k) )

          is_on_same_level =                                            &
                  (hhl(i  ,j  ,k+1)<=zu) .AND. (zu <= hhl(i  ,j  ,k) )  &  
            .AND. (hhl(i+1,j  ,k+1)<=zu) .AND. (zu <= hhl(i+1,j  ,k) )  &  
            .AND. (hhl(i  ,j  ,k+1)<=zv) .AND. (zv <= hhl(i  ,j  ,k) )  &  
            .AND. (hhl(i  ,j+1,k+1)<=zv) .AND. (zv <= hhl(i  ,j+1,k) )  

          IF ( .NOT. is_on_same_level ) THEN
            is_search_finished = .TRUE.  
            EXIT               
          END IF

        END DO
        IF ( is_search_finished ) EXIT
      END DO
      IF ( is_search_finished ) EXIT

    END DO
    k_index_same_level_Mahrer = k-1

    ! At least necessary for domain decomposition independency
    ! (or processor reproducibility):
    CALL global_values( k_index_same_level_Mahrer, 1, 'MIN', imp_integers, &
      &                 icomm_cart ,-1, yzerrmsg, izerror )

    DEALLOCATE( hml,    STAT=istat )
    DEALLOCATE( k_star, STAT=istat  )

  END IF

  ! lower boundary condition for w
  !----------------------------------------------------------------------

  ! a.) extrapolation order of u, v to the bottom:

  IF (      (itype_bbc_w==1) .OR.     &
            (itype_bbc_w==3) .OR.     &
            (itype_bbc_w==5) .OR.     &
            (itype_bbc_w==7) .OR.     &
            (itype_bbc_w==9) ) THEN

    ! free slip condition (no extrapolation):
    !i_bbc_w_extrapol_uv_order = 0

    ! linear extrapolation
    !    (levels ke-1 and ke to ke+1/2)
    i_bbc_w_extrapol_uv_order = 1

  ELSE IF ( (itype_bbc_w==0) .OR.     &
            (itype_bbc_w==2) .OR.     &
            (itype_bbc_w==4) .OR.     &
            (itype_bbc_w==6) .OR.     &
            (itype_bbc_w==8) ) THEN

    ! quadratic extrapolation (by G. Zaengl)
    !    (levels ke-2, ke-1 and ke to ke+1/2)
    i_bbc_w_extrapol_uv_order = 2

  ELSE IF ( itype_bbc_w >= 100 ) THEN
    ! alternative nomenclature for itype_bbc_w:
    ! "1ed" : e = extrapolation order of u, v to the bottom surface (0,1,or 2)
    !         d = discretization order of dh/dx, dh/dy (currently 2,3,or 4)

    IF (      MOD( itype_bbc_w / 10, 10) == 0 ) THEN
      ! free slip condition (no extrapolation):
      i_bbc_w_extrapol_uv_order = 0

    ELSE IF ( MOD( itype_bbc_w / 10, 10) == 1 ) THEN
      ! linear extrapolation:
      i_bbc_w_extrapol_uv_order = 1

    ELSE IF ( MOD( itype_bbc_w / 10, 10) == 2 ) THEN
      ! quadratic extrapolation:
      i_bbc_w_extrapol_uv_order = 2

    ELSE
      yzerrmsg = "false value in itype_bbc_w"
      CALL model_abort (my_cart_id, 1234, yzerrmsg, "init_fast_waves_sc_1")
    END IF

  ELSE
    yzerrmsg = "false value in itype_bbc_w"
    CALL model_abort (my_cart_id, 1234, yzerrmsg, "init_fast_waves_sc_1")
  END IF

  ! b.) discretization order of dh/dx and dh/dy:

  IF ( itype_bbc_w >= 100 ) THEN
    i_bbc_w_dhdx_order = MOD( itype_bbc_w, 10)

    IF ( ( i_bbc_w_dhdx_order < 2 ) .OR. ( i_bbc_w_dhdx_order > 4 ) ) THEN
      yzerrmsg = "false value in i_bbc_w_dhdx_order: is not equal 2,3, or 4!"
      CALL model_abort (my_cart_id, 1234, yzerrmsg, "init_fast_waves_sc_1")
    END IF

  ELSE 
    ! old nomenclature for itype_bbc_w is used; proceed with 
    i_bbc_w_dhdx_order = -1  ! = 'undefined'
  END IF


  IF ( (itype_bbc_w >= 100) .AND. (itype_bbc_w /= 114) ) THEN
    IF ( my_cart_id == 0 ) THEN
      WRITE(*,*) "WARNING: you don't use the recommended value itype_bbc_w=114!"
    END IF
  END IF

  ! Other metric terms:
  ! (in the case of memory problems during physics calculations
  ! the following part might be put into init_fast_waves_sc_2)
  !-----------------------------------------------------------------

  ALLOCATE( dz_dlam  (ie, je, ke1), STAT=istat)
  ALLOCATE( dz_dphi  (ie, je, ke1), STAT=istat)

  ALLOCATE( dz_lam_u_half (ie, je, ke), STAT=istat)
  ALLOCATE( dz_phi_v_half (ie, je, ke), STAT=istat)


  IF (     (itype_bbc_w == 4) .OR. (itype_bbc_w == 5)    &
      .OR. (i_bbc_w_dhdx_order == 4 ) ) THEN

    ! additionally calculate
    ! 4th order approximation of dz/dlam, dz/dphi
    ! (only used for the free slip bottom BC for w)

    ALLOCATE( dz_dlam_4ord (ie,je), STAT=istat )
    ALLOCATE( dz_dphi_4ord (ie,je), STAT=istat )

    CALL calc_dz_dx( dz_dlam, dz_dphi,               &
                     dz_lam_u_half, dz_phi_v_half,   &
                     dz_dlam_4ord, dz_dphi_4ord )

  ELSE

    CALL calc_dz_dx( dz_dlam, dz_dphi,               &
                     dz_lam_u_half, dz_phi_v_half )

  END IF

END SUBROUTINE init_fast_waves_sc_1

!=============================================================================
!=============================================================================

SUBROUTINE init_fast_waves_sc_2

  !---------------------------------------------------------------------------
  !
  ! Description:
  !   Initialization subroutine (2).
  !   It is called at every LARGE timestep.
  !
  !---------------------------------------------------------------------------

  INTEGER :: i,j,k

  DO k=1, ke
    DO j=1, je
      DO i=1, ie
        sqrtgtot_r(i,j,k) = acrlat(j,1) * sqrtg_r_s(i,j,k)
      END DO
    END DO
  END DO

END SUBROUTINE init_fast_waves_sc_2

!=============================================================================
!=============================================================================

SUBROUTINE init_fast_waves_sc_3( dtsmall )

  !--------------------------------------------------------------------------
  !
  ! Description:
  !   Initialization subroutine (3).
  !   Like "init_fast_waves_sc_1", it is called only once at the beginning of 
  !   the simulation, but here the small time step must be known.
  !
  !--------------------------------------------------------------------------

  REAL (KIND=ireals), INTENT(IN) :: dtsmall

  ! Initialization of a possible Rayleigh damping by Klemp et al. (2008)
  !----------------------------------------------------------------------

  INTEGER :: k

  ALLOCATE(wrdfac(ke))

  ! Damping factor for new sponge condition (Klemp et al., 2008, MWR)
  IF (lspubc .AND. (itype_spubc == 3) ) THEN
    DO k = 1, ke_rd
      wrdfac(k) = 1.0_ireals/(1.0_ireals + 10.0_ireals*rdcoef(k)*dtsmall)
    ENDDO
  ENDIF

END SUBROUTINE init_fast_waves_sc_3

!=============================================================================
!=============================================================================

SUBROUTINE finalize_fast_waves_sc

  !---------------------------------------------------------------------------
  !
  ! Description:
  !   proper finalization of 'fast_waves_sc'
  !   (at the end of the whole simulation).
  !
  !---------------------------------------------------------------------------

  IF ( l_hor_p_grad_Mahrer ) THEN
    DEALLOCATE( p_hor_idx )
    DEALLOCATE( p_hor_wght )
    DEALLOCATE( k_star_u )
    DEALLOCATE( k_star_v )
  END IF

  DEALLOCATE( wrdfac )

  DEALLOCATE( dz_dlam )
  DEALLOCATE( dz_dphi )

  DEALLOCATE( dz_lam_u_half )
  DEALLOCATE( dz_phi_v_half )

  IF (     (itype_bbc_w == 4) .OR. (itype_bbc_w == 5)    &
      .OR. (i_bbc_w_dhdx_order == 4 ) ) THEN

    DEALLOCATE( dz_dlam_4ord )
    DEALLOCATE( dz_dphi_4ord )
  END IF

END SUBROUTINE finalize_fast_waves_sc

!=============================================================================
!=============================================================================

SUBROUTINE alloc_fast_waves_sc

  !---------------------------------------------------------------------------
  !
  ! Description:
  !   Allocation of variables OUTSIDE of fast_waves_sc
  !
  !---------------------------------------------------------------------------
  INTEGER :: izstata

  ALLOCATE( lgs_A         (ie, je, ke1), STAT=izstata )
  ALLOCATE( lgs_B         (ie, je, ke1), STAT=izstata )
  ALLOCATE( lgs_C         (ie, je, ke1), STAT=izstata )

  ALLOCATE( T_tot         (ie, je, ke ), STAT=izstata )
  ALLOCATE( gamma_p_tot   (ie, je, ke ), STAT=izstata )  ! = c_p/c_v*p_tot
  ALLOCATE( rho           (ie, je, ke ), STAT=izstata )
  ALLOCATE( rho_inv       (ie, je, ke ), STAT=izstata )
  ALLOCATE( rho_inv_k_avg (ie, je, ke ), STAT=izstata )

  ALLOCATE( p0T_pT0_k_avg (ie, je, 2:ke), STAT=izstata ) ! = (p0*T)/(p*T0)
  ALLOCATE( p0_pT0_k_avg  (ie, je, 2:ke), STAT=izstata ) ! = p0/(p*T0)
  ALLOCATE( p_inv_k_avg   (ie, je, 2:ke), STAT=izstata ) ! = 1/p

  ALLOCATE( frac_1_rhodx  (ie, je, ke ), STAT=izstata )
  ALLOCATE( frac_1_rhody  (ie, je, ke ), STAT=izstata )

  ALLOCATE( sqrtgtot_r    (ie, je, ke ), STAT=izstata)

  ALLOCATE( alpha_div_h   (ie, je, ke ), STAT=izstata )

END SUBROUTINE alloc_fast_waves_sc

!=============================================================================
!=============================================================================

SUBROUTINE dealloc_fast_waves_sc

  !---------------------------------------------------------------------------
  !
  ! Description:
  !   Deallocation of variables OUTSIDE of fast_waves_sc
  !
  !---------------------------------------------------------------------------

  DEALLOCATE( lgs_A )
  DEALLOCATE( lgs_B )
  DEALLOCATE( lgs_C )

  DEALLOCATE( T_tot )
  DEALLOCATE( gamma_p_tot   )
  DEALLOCATE( rho     )
  DEALLOCATE( rho_inv )
  DEALLOCATE( rho_inv_k_avg )

  DEALLOCATE( p0T_pT0_k_avg )
  DEALLOCATE( p0_pT0_k_avg )
  DEALLOCATE( p_inv_k_avg )

  DEALLOCATE( frac_1_rhodx )
  DEALLOCATE( frac_1_rhody )

  DEALLOCATE( sqrtgtot_r )

  DEALLOCATE( alpha_div_h )

END SUBROUTINE dealloc_fast_waves_sc

!=============================================================================
!=============================================================================

SUBROUTINE fast_waves_strong_conserv ( n_small_steps, dts, xkd,             &
  u, v, w, Tp, pp, qv, q_cond,                                  &
  du_dt_slow, dv_dt_slow, dw_dt_slow, dtp_dt_slow, dpp_dt_slow, &
  zubdt_west, zubdt_east, zvbdt_north, zvbdt_south )

  !---------------------------------------------------------------------------
  !
  ! Description:
  !   integrate the sound expansion and buoyancy/gravity wave expansion terms
  !   over n_small_steps small time steps dts.
  !   In every such small time step integration the (constant) slow tendencies
  !   are added.
  !   
  !   This module procedure calculates the new values of the prognostic 
  !   variables  u, v, w, pp and Tp.
  !
  ! input:
  !   n_small_steps : number of small time steps
  !   dts           : small time step (in RK: dts = dt / n_steps )
  !   xkd           : dim.-less coefficient for divergence damping
  !   u, v, w, Tp, pp, qv, q_cond  : fields of the 'dynamical core'
  !   du_dt_slow, dv_dt_slow, ...  : slow tendencies
  !
  ! output:
  !   u, v, w, Tp, pp: updated fields 
  !
  ! Method:
  !   - horizontally forward-backward, vertically implicit (Crank-Nicholson)
  !     perform the solution of the vertically implicit equation for w
  !       A(k) w(k-1) + B(k) w(k) + C(k) w(k+1) = rhs(k)
  !              ( half levels w[k] -> w(zeta = k-1/2) )
  !   - correct weighting of all vertical derivatives and averages
  !   - strong conservation form of the divergence operator
  !   - optional: use of isotropic divergence damping
  !   - optional: Mahrer (1984) discret. of horizontal pressure gradients
  !
  !---------------------------------------------------------------------------

  ! Subroutine arguments:
  ! --------------------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
    n_small_steps         ! number of small time steps for integration

  REAL    (KIND=ireals   ), INTENT (IN) ::  &
    dts,                & ! small time step
    xkd                   ! constant factor for divergence-type damping

  REAL    (KIND=ireals   ), INTENT (INOUT) ::  &
    u (ie, je, ke),     & !
    v (ie, je, ke),     & !
    w (ie, je, ke1),    & !
    Tp(ie, je, ke),     & !
    pp(ie, je, ke)        !

  REAL    (KIND=ireals   ), INTENT (IN) ::  &
    qv     (ie, je, ke), & !
    q_cond (ie, je, ke)    !

  REAL    (KIND=ireals   ), INTENT (IN) ::  &
    du_dt_slow  (ie,je,ke),   & ! slow tendency of u-velocity
    dv_dt_slow  (ie,je,ke),   & ! slow tendency of v-velocity
    dw_dt_slow  (ie,je,ke1),  & ! slow tendency of w-velocity
    dtp_dt_slow (ie,je,ke),   & ! slow tendency of temperature
    dpp_dt_slow (ie,je,ke)      ! slow tendency of pressure perturbation

  REAL    (KIND=ireals   ), INTENT (IN) ::  &
    zubdt_west (je,ke), & ! u-boundary tendency at west boundary (total domain)
    zubdt_east (je,ke), & ! u-boundary tendency at east boundary (total domain)
    zvbdt_south(ie,ke), & ! v-boundary tendency at south boundary (total domain)
    zvbdt_north(ie,ke)    ! v-boundary tendency at north boundary (total domain)


  ! Local scalars:
  ! -----------------------------
  INTEGER (KIND=iintegers) ::  &
    kzdims(24),          & ! Vertical dimensions for exchg_boundaries
    i,  j,  k,           & ! Loop indices
    k_top, k_bottom,     &
    izerror,             & ! error status
    izstata,             & ! error status at allocation
    izstatd                ! error status at deallocation

  CHARACTER (LEN=80)        :: yzerrmsg


  ! Local (allocatable) arrays:
  ! -----------------------------

  ! REAL(KIND=ireals), ALLOCATABLE ::  zhml(:,:,:)

  ! automatic fields are a bit faster than allocations:
  REAL (KIND=ireals) :: Z_horiz            (ie, je, ke1)
  REAL (KIND=ireals) :: div_tilde_horiz    (ie, je, ke)
  REAL (KIND=ireals) :: div_tilde_vertic_uv(ie, je, ke)

  REAL (KIND=ireals) :: p0_inv             (ie, je, ke)

  !REAL (KIND=ireals) :: div_expl_p         (ie, je, ke)
  !REAL (KIND=ireals) :: div_expl_T         (ie, je, ke)
  REAL (KIND=ireals) :: div_expl_p_T       (ie, je, ke)  ! beta_s_6=beta_s_4 necessary
  REAL (KIND=ireals) :: div_expl_uv        (ie, je, ke)
  REAL (KIND=ireals), ALLOCATABLE :: div_expl_w (:,:,:)

  REAL (KIND=ireals) :: pp_minus_div       (ie, je, ke)
  REAL (KIND=ireals) :: w_k_avg_old        (ie, je, ke)
  REAL (KIND=ireals) :: q_x_k_avg          (ie, je, 2:ke)
  REAL (KIND=ireals) :: rho0_rho_q_x_k_avg (ie, je, 2:ke)

  REAL (KIND=ireals) :: b_4p               (ie, je, ke)

  REAL (KIND=ireals) :: sqrtg_eddlon_u       (ie, je, ke)
  REAL (KIND=ireals) :: sqrtg_cosphi_eddlat_v(ie, je, ke)

  REAL (KIND=ireals) :: contrib_u          (ie, je)
  REAL (KIND=ireals) :: contrib_v          (ie, je)

  ! r_earth * cos(phi) (at scalar position):
  REAL (KIND=ireals) :: r_cosphi           (je)

  ! other auxiliary fields:
  REAL (KIND=ireals) :: work3              (ie, ke)

  REAL (KIND=ireals), ALLOCATABLE :: dpp_dz_s (:,:,:)

  REAL (KIND=ireals), ALLOCATABLE :: u_sfc(:,:), v_sfc(:,:)

  REAL (KIND=ireals) :: buoy
  REAL (KIND=ireals) :: dug_dlam
  REAL (KIND=ireals) :: dvg_dphi
  REAL (KIND=ireals) :: dw_dzeta
  REAL (KIND=ireals) :: div_tilde_vertic_w, div_tilde_vertic
  !REAL (KIND=ireals) :: div_part_p, div_part_T  ! use, if beta_s_9 /= beta_s_8
  REAL (KIND=ireals) :: div_part_p_T             ! beta_s_9=beta_s_8 necessary
  REAL (KIND=ireals) :: div_impl_T 
  REAL (KIND=ireals) :: dpp_dlam_metr, dpp_dphi_metr
  REAL (KIND=ireals) :: Tp_k_avg
  REAL (KIND=ireals) :: pp_k_avg
  REAL (KIND=ireals) :: w_k_avg_new
  REAL (KIND=ireals) :: contr_b_4p, contr_b_T
  REAL (KIND=ireals) :: p_eq_expl, T_eq_expl, b_w
  REAL (KIND=ireals) :: one_m_wgtfac
  REAL (KIND=ireals) :: dpp_dz_u_double, dpp_dz_v_double
  REAL (KIND=ireals) :: eddlat_cosphi

  REAL (KIND=ireals) :: r_earth_recip
  REAL (KIND=ireals) :: dlambda, dphi   ! in radians
  REAL (KIND=ireals) :: dts_inv

  ! auxiliary variables for the Mahrer (1984)/Zaengl (2012) discretization
  ! of horizontal pressure gradients
  REAL (KIND=ireals) :: buoy_r_o, buoy_r_u, buoy_l_o, buoy_l_u
  REAL (KIND=ireals) :: buoy_r, buoy_l
  REAL (KIND=ireals) :: p_r, p_l
  REAL (KIND=ireals) :: zu_k_star, zv_k_star
  REAL (KIND=ireals), ALLOCATABLE :: dp_dlam_k_star(:,:)
  REAL (KIND=ireals), ALLOCATABLE :: dp_dphi_k_star(:,:)
  REAL (KIND=ireals), ALLOCATABLE :: dbuoy_dlam_k_star(:,:)
  REAL (KIND=ireals), ALLOCATABLE :: dbuoy_dphi_k_star(:,:)

  REAL (KIND=ireals) :: zu, zv
  REAL (KIND=ireals) :: z1, z2

  INTEGER :: n_step

  REAL (KIND=ireals) :: switch_div_damp

  INTEGER :: istart_1, iend_1
  INTEGER :: jstart_1, jend_1

  INTEGER :: m
  INTEGER :: k1, k2

  ! boundary fields for the radiation condition:
  REAL (KIND=ireals), ALLOCATABLE ::  bd_west  (:,:,:)
  REAL (KIND=ireals), ALLOCATABLE ::  bd_east  (:,:,:)
  REAL (KIND=ireals), ALLOCATABLE ::  bd_north (:,:,:)
  REAL (KIND=ireals), ALLOCATABLE ::  bd_south (:,:,:)

  REAL (KIND=ireals) :: dp_dx_rho, dp_dy_rho

  !
  ! End of header
  !===========================================================================

  !---------------------------------------------------------------------------
  ! Begin Subroutine fast_waves_strong_conserv
  !---------------------------------------------------------------------------

  ! Initialize, whether debug output shall be done
  ! ----------------------------------------------------------------

  IF (ldebug_dyn) THEN
    IF (lprintdeb_all) THEN
      izdebug = idbg_level
    ELSE
      IF (my_cart_id == 0) THEN
        izdebug = idbg_level
      ELSE
        izdebug = 0
      ENDIF
    ENDIF
  ELSE
    izdebug = 0
  ENDIF

  IF (izdebug >= 10) THEN
    WRITE(*,*) "Subr. [fast_waves_strong_conserv]..."
  END IF

  ! Allocations
  !-----------------------------------------------------------------

  ! (allocations instead of automatic fields are used, because either 
  !  fields are declared outside of this subroutine 
  !  or use of those fields is case dependent)


  ALLOCATE( lgs_rhs  (ie, je, ke1), STAT=izstata)

  ALLOCATE( u_sfc    (ie,je), STAT=izstata)
  ALLOCATE( v_sfc    (ie,je), STAT=izstata)

  IF ( .NOT. l_calc_lhs_at_1st_RKstep ) THEN
    CALL alloc_fast_waves_sc 
    CALL init_fast_waves_sc_2
  END IF

  IF ( l_3D_div_damping ) THEN
    ALLOCATE( div_expl_w (ie, je, ke),  STAT=izstata)
  END IF

  IF ( l_hor_p_grad_Mahrer ) THEN

    ALLOCATE( dp_dlam_k_star    (ie, je), STAT=izstata)
    ALLOCATE( dp_dphi_k_star    (ie, je), STAT=izstata)
    ALLOCATE( dbuoy_dlam_k_star (ie, je), STAT=izstata)
    ALLOCATE( dbuoy_dphi_k_star (ie, je), STAT=izstata)

  END IF

  ALLOCATE( dpp_dz_s(ie, je, vcoord%kflat:ke), STAT=izstata)

  IF ( lradlbc ) THEN
    ALLOCATE (  bd_west ( je, ke1, 1:5 ) ) ! 1=u, 2=v, 3=w, 4=Tp, 5=pp
    ALLOCATE (  bd_east ( je, ke1, 1:5 ) )
    ALLOCATE (  bd_north( ie, ke1, 1:5 ) )
    ALLOCATE (  bd_south( ie, ke1, 1:5 ) )
  END IF

  ! --- set index bounds for u and v integration ---

  ! western boundary:
  IF ( my_cart_neigh(1) == -1 ) THEN
    i_start_u = istartu
  ELSE
    i_start_u = istartu-1
  END IF

  ! southern boundary:
  IF ( my_cart_neigh(4) == -1 ) THEN
    j_start_v = jstartv
  ELSE
    j_start_v = jstartv-1
  END IF

  istart_1 = MIN( i_start_u, istartv   )
  iend_1   = MAX( iendu,     iendv     )

  jstart_1 = MIN( jstartu,   j_start_v )
  jend_1   = MAX( jendu,     jendv     )


  ! some abbreviations for runtime optimizations
  !-----------------------------------------------------------------

  r_earth_recip = 1.0_ireals / r_earth
  dts_inv       = 1.0_ireals / dts

  p0_inv(:,:,:) = 1.0_ireals / p0(:,:,:)

  DO j=1, je
    r_cosphi(j) = r_earth * crlat(j,1)
  END DO

  DO k=1, ke
    DO j=1, je
      eddlat_cosphi = eddlat * crlat(j,2)
      DO i=1, ie
        sqrtg_eddlon_u       (i,j,k) = eddlon        / sqrtg_r_u(i,j,k)
        sqrtg_cosphi_eddlat_v(i,j,k) = eddlat_cosphi / sqrtg_r_v(i,j,k)
      END DO
    END DO
  END DO

  DO j=jstart, jend

    !CDIR ON_ADB(work3)
    DO k=1, ke
      DO i=istart, iend
        work3(i,k) = rvd_m_o * qv(i,j,k) - q_cond(i,j,k)
      END DO
    END DO

    DO k=2, ke
      DO i=istart, iend
        q_x_k_avg(i,j,k) =       wgtfac(i,j,k)  * work3(i,k)    &
          &        + (1.0_ireals-wgtfac(i,j,k)) * work3(i,k-1)
      END DO
    END DO

  END DO

  ! divergence damping coefficients:
  IF ( .NOT. l_calc_lhs_at_1st_RKstep ) THEN
    CALL init_div_damping_coeff( xkd, dts, alpha_div_h, alpha_div_v_to_h )
  END IF

  IF ( (.NOT. l_calc_lhs_at_1st_RKstep) .AND. l_small_pert_in_pT ) THEN
    ! assumption: p and T in all coefficients change only very little 
    ! during the small time step integration
    CALL calc_total_T_p_rho( Tp, pp, qv, q_cond)
  END IF

  IF ( l_small_pert_in_pT ) THEN
    rho0_rho_q_x_k_avg (istart:iend, jstart:jend, :) =   & 
      &   p0T_pT0_k_avg(istart:iend, jstart:jend, :)     &
      &     * q_x_k_avg(istart:iend, jstart:jend, :) 

  END IF

  div_tilde_vertic_uv(:,:,1:MAX(1,vcoord%kflat-1)) = 0.0_ireals

  div_tilde_horiz(:,:,1)  = 0.0_ireals
  div_tilde_horiz(:,:,ke) = 0.0_ireals

  ! set boundary, to avoid access to NaN (??):
  Z_horiz(1:istart,:,:) = 0.0_ireals
  Z_horiz( iend:ie,:,:) = 0.0_ireals
  Z_horiz(:,1:jstart,:) = 0.0_ireals
  Z_horiz(:, jend:je,:) = 0.0_ireals

  IF ( .NOT. l_hor_p_grad_Mahrer ) THEN
    dpp_dz_s(1:istart+2,:,:) = 0.0_ireals
    dpp_dz_s(iend-2:ie ,:,:) = 0.0_ireals
    dpp_dz_s(:,1:jstart+2,:) = 0.0_ireals
    dpp_dz_s(:,jend-2:je ,:) = 0.0_ireals
  END IF

  div_expl_uv (:,:,:)     = 0.0_ireals

  !---------------------------------------------------------------------------
  !   small time step integration
  !---------------------------------------------------------------------------

  DO n_step = 1, n_small_steps 

    IF ( .NOT. l_small_pert_in_pT ) THEN
      ! assumption: the change of p and T in the coefficients is significant
      ! during the small time step integration
      CALL calc_total_T_p_rho( Tp, pp, qv, q_cond )

      rho0_rho_q_x_k_avg (istart:iend, jstart:jend, :) =    &  
        &   p0T_pT0_k_avg(istart:iend, jstart:jend, :)      &
        &     * q_x_k_avg(istart:iend, jstart:jend, :) 

    END IF


    IF ( l_use_lateral_bound_values ) THEN
      ! Set boundary values for u and v from the boundary tendencies
      ! to compute the wind divergence at the scalar gridpoints along
      ! the boundary lines

      IF ( .NOT.lradlbc .AND. .NOT.lperi_x) THEN
        ! west
        IF (my_cart_neigh(1) == -1) THEN
          DO  k = 1, ke
            !CDIR NOLOOPCHG
            DO  j = jstartu, jendu
              u(istartu-1,j,k) = u(istartu-1,j,k) + zubdt_west(j,k)
            ENDDO
          END DO
        ENDIF
        ! east
        IF (my_cart_neigh(3) == -1) THEN
          DO  k = 1, ke
            !CDIR NOLOOPCHG
            DO  j = jstartu, jendu
              u(iendu+1,j,k) = u(iendu+1,j,k) + zubdt_east(j,k)
            ENDDO
          END DO
        ENDIF
      END IF

      IF ( .NOT.lradlbc  .AND. .NOT.lperi_y) THEN
        ! south
        IF (my_cart_neigh(4) == -1) THEN
          DO  k = 1, ke
            !CDIR NOLOOPCHG
            DO i = istartv, iendv
              v(i,jstartv-1,k) = v(i,jstartv-1,k) + zvbdt_south(i,k)
            ENDDO
          END DO
        ENDIF
        ! north
        IF (my_cart_neigh(2) == -1) THEN
          DO  k = 1, ke
            !CDIR NOLOOPCHG
            DO i = istartv, iendv
              v(i,jendv+1,k) = v(i,jendv+1,k) + zvbdt_north(i,k)
            ENDDO
          END DO
        ENDIF
      END IF

    ELSE

      IF ( lradlbc .AND. .NOT.( lperi_x .OR. lperi_y ) ) THEN
        ! radiative lateral boundary condition (init)

        CALL radiative_lbc(u, v, w, pp, Tp, dts,                    &
                           bd_west, bd_east, bd_north, bd_south )
        ! remark: the update of the boundary fields with bd_west, ...
        ! takes place at the end of the time loop

      END IF

    END IF


    ! calculate the explicit parts of the divergence
    !---------------------------------------------------------------

    IF ( n_step == 1 ) THEN
      ! these calculations take place also later on for u^(n+1) and v^(n+1):

      CALL calc_uv_surface( u, v, i_bbc_w_extrapol_uv_order, u_sfc, v_sfc)
      CALL calc_Z_horiz( u, v, u_sfc, v_sfc, Z_horiz )

      DO k=vcoord%kflat, ke
        DO j=jstart, jend
          DO i=istart, iend
            div_tilde_vertic_uv(i,j,k) = Z_horiz(i,j,k+1) - Z_horiz(i,j,k)
          END DO
        END DO
      END DO

      IF ( l2dim ) THEN
        DO j=1, nboundlines
          div_tilde_vertic_uv(:,jstart-j,:) = div_tilde_vertic_uv(:,jstart,:)
          div_tilde_vertic_uv(:,jend  +j,:) = div_tilde_vertic_uv(:,jstart,:)
        END DO
      END IF


      IF ( .NOT. l_no_div_damping_in_1st_step ) THEN

        DO k=1, ke
          DO j=jstart, jend
            DO i=istart, iend

              dug_dlam  = u(i  ,j,k) * sqrtg_eddlon_u(i,  j,k)           &
                &       - u(i-1,j,k) * sqrtg_eddlon_u(i-1,j,k)

              dvg_dphi  = v(i,j  ,k) * sqrtg_cosphi_eddlat_v(i,j  ,k)    &  
                &       - v(i,j-1,k) * sqrtg_cosphi_eddlat_v(i,j-1,k)

              div_tilde_horiz(i,j,k) = dug_dlam + dvg_dphi

            END DO
          END DO
        END DO

        IF ( l2dim ) THEN
          DO j=1, nboundlines
            div_tilde_horiz(:,jstart-j,:) = div_tilde_horiz(:,jstart,:)
            div_tilde_horiz(:,jend  +j,:) = div_tilde_horiz(:,jstart,:)
          END DO
        END IF

      END IF

    END IF     ! n_step==1

    ! in principle: the following could be used also for k=1 and k=ke.
    ! But one can use the BC dzeta/dt=0 if:
    !    beta_s_9=beta_s_6
    !    beta_s_8=beta_s_4
    !    beta_d_6=beta_d_4

    IF ( i_div_at_horiz_bound == 1 ) THEN
      ! extra consideration for k=1 and k=ke
      ! e.g. to use d zeta/dt=0 at 'horizontal' boundaries:
      k_top    = 2
      k_bottom = ke-1
    ELSE IF ( i_div_at_horiz_bound == 2 ) THEN
      k_top    = 1
      k_bottom = ke
    END IF

    IF ( ( n_step==1 ) .AND. l_no_div_damping_in_1st_step ) THEN

      DO k=k_top, k_bottom
        DO j=jstart, jend
          DO i=istart, iend

            div_tilde_vertic_w  = -r_cosphi(j) * ( w(i,j,k+1) - w(i,j,k) )

            div_expl_p_T(i,j,k) = sqrtgtot_r(i,j,k) * (        &
                  !beta_s_3_expl * div_tilde_horiz(i,j,k)      & ! beta_s_3_expl=0
              &    beta_s_8_expl * div_tilde_vertic_uv(i,j,k)  &
              &  + beta_s_4_expl * div_tilde_vertic_w ) 

          END DO
        END DO
      END DO

    ELSE   !  IF  ( n_step>1 ) .OR. (.NOT. l_no_div_damping_in_1st_step)

      DO k=k_top, k_bottom

        IF ( k < vcoord%kflat ) THEN

          DO j=jstart, jend
            DO i=istart, iend

              div_tilde_vertic_w  = -r_cosphi(j) * ( w(i,j,k+1) - w(i,j,k) )

              div_expl_p_T(i,j,k) = sqrtgtot_r(i,j,k) * (   &
                    !beta_s_3_expl * div_tilde_horiz(i,j,k) & !beta_s_3_expl=0
                &  + beta_s_4_expl * div_tilde_vertic_w ) 

              div_expl_uv(i,j,k) = sqrtgtot_r(i,j,k) * ( &
                &    div_tilde_horiz(i,j,k)              & !*beta_d_1_expl(=1)
                &  + div_tilde_vertic_w )                  !*beta_d_2_expl(=1)

            END DO
          END DO

        ELSE   ! case k >= vcoord%kflat

          DO j=jstart, jend
            DO i=istart, iend

              div_tilde_vertic_w  = -r_cosphi(j) * ( w(i,j,k+1) - w(i,j,k) )

              div_expl_p_T(i,j,k) = sqrtgtot_r(i,j,k) * (        &
                    !beta_s_3_expl * div_tilde_horiz(i,j,k)      & ! beta_s_3_expl=0
                &    beta_s_8_expl * div_tilde_vertic_uv(i,j,k)  &
                &  + beta_s_4_expl * div_tilde_vertic_w ) 

              div_expl_uv(i,j,k) = sqrtgtot_r(i,j,k) * (  &
                &    div_tilde_horiz(i,j,k)               & ! * beta_d_1_expl (=1)
                &  + div_tilde_vertic_uv(i,j,k)           & ! * beta_d_5_expl (=1)
                &  + div_tilde_vertic_w )                   ! * beta_d_2_expl (=1)

            END DO
          END DO

        END IF

      END DO

    END IF

    IF ( l_3D_div_damping .AND. &
           (      (       l_no_div_damping_in_1st_step .AND. ( n_step>1 ) )   &
             .OR. ( .NOT. l_no_div_damping_in_1st_step )    )                 &
       ) THEN
      DO k=k_top, k_bottom
        DO j=jstart, jend
          DO i=istart, iend

            div_tilde_vertic_w  = -r_cosphi(j) * ( w(i,j,k+1) - w(i,j,k) )

            div_expl_w (i,j,k) = sqrtgtot_r(i,j,k) * (         &
              &    beta_d_3_expl * div_tilde_horiz(i,j,k)      &
              &  + beta_d_6_expl * div_tilde_vertic_uv(i,j,k)  &
              &  + beta_d_4_expl * div_tilde_vertic_w ) 

          END DO
        END DO
      END DO
    END IF


    IF ( i_div_at_horiz_bound == 1 ) THEN

      ! top layer:

      k=1

      DO j=jstart, jend
        DO i=istart, iend

          div_tilde_vertic = Z_horiz(i,j,2) - r_cosphi(j) * w(i,j,2)

          div_expl_p_T(i,j,k) = sqrtgtot_r(i,j,k) * (  &
              !beta_s_3_expl * div_tilde_horiz(i,j,k)  & ! beta_s_3_expl=0
               beta_s_4_expl * div_tilde_vertic )        ! beta_s_8=beta_s_4 necessary 

          div_expl_uv(i,j,k) = sqrtgtot_r(i,j,k) * (   &
                               div_tilde_horiz(i,j,k)  & ! * beta_d_1_expl (=1)
                             + div_tilde_vertic )        ! beta_d_5=beta_d_2 necessary

        END DO
      END DO

      IF ( l_3D_div_damping .AND.                                              &
             (      (       l_no_div_damping_in_1st_step .AND. ( n_step>1 ) )  &
               .OR. ( .NOT. l_no_div_damping_in_1st_step )    )                &
         ) THEN

        DO j=jstart, jend
          DO i=istart, iend

            ! div_tilde_vertic = (r*cos*g*wcon(i,j,2) - r*cos*g*wcon(i,j,1) )/dzeta:
            !                                             =0!
            div_tilde_vertic = Z_horiz(i,j,2) - r_cosphi(j) * w(i,j,2)

            div_expl_w (i,j,k) = sqrtgtot_r(i,j,k) * (      &
                  beta_d_3_expl * div_tilde_horiz(i,j,k)    &
                + beta_d_4_expl * div_tilde_vertic )     ! beta_d_6=beta_d_4 necessary  
          END DO
        END DO
      END IF


      ! bottom layer:

      k=ke

      DO j=jstart, jend
        DO i=istart, iend

          ! div_tilde_vertic = [ Z(i,j,ke+1) - Z(i,j,ke) ] / dzeta:
          !                         =0!                       =1

          div_tilde_vertic = - ( Z_horiz(i,j,ke) - r_cosphi(j) * w(i,j,ke) )

          div_expl_p_T(i,j,k) = sqrtgtot_r(i,j,k) * (   &
              !beta_s_3_expl * div_tilde_horiz(i,j,k)   & ! beta_s_3_expl=0
               beta_s_4_expl * div_tilde_vertic )         ! beta_s_8=beta_s_4 necessary 

          div_expl_uv(i,j,k) = sqrtgtot_r(i,j,k) * (   &
                               div_tilde_horiz(i,j,k)  & ! * beta_d_1_expl (=1)
                             + div_tilde_vertic )        ! beta_d_5=beta_d_2 necessary

        END DO
      END DO

      IF ( l_3D_div_damping .AND.                                              &
             (      (       l_no_div_damping_in_1st_step .AND. ( n_step>1 ) )  &
               .OR. ( .NOT. l_no_div_damping_in_1st_step )    )                &
         ) THEN

        DO j=jstart, jend
          DO i=istart, iend

            ! div_tilde_vertic= [ Z(i,j,ke+1) - Z(i,j,ke) ] / dzeta :
            !                         =0!                      =1
            div_tilde_vertic = - ( Z_horiz(i,j,ke) - r_cosphi(j) * w(i,j,ke) )

            div_expl_w (i,j,k) = sqrtgtot_r(i,j,k) * (      &
                   beta_d_3_expl * div_tilde_horiz(i,j,k)   &
                 + beta_d_4_expl * div_tilde_vertic )    ! beta_d_6=beta_d_4 necessary

          END DO
        END DO

      END IF

    END IF

    ! ----- Exchange of pp and div_expl_uv: ---------------

    kzdims(1:24) =                        &
         (/  ke, ke ,   0,   0,   0,      &
              0,   0,   0,   0,   0,      &
              0,   0,   0,   0,   0,      &
              0,   0,   0,   0,   0,      &
              0,   0,   0,   0  /)

    ! non-optimized version (icase=0 <--> ldatatypes=.FALSE.)
    CALL exchg_boundaries                                                    &
        (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie,je, &
         kzdims, jstartpar, jendpar, nboundlines, nboundlines,my_cart_neigh, &
         lperi_x, lperi_y, l2dim,                                            &
         10000+nexch_tag, .FALSE., ncomm_type, izerror, yzerrmsg,            &
         pp         (:,:,:),           &
         div_expl_uv(:,:,:)  )

    IF ( izerror /= 0_iintegers ) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'fast_waves_strong_conserv')
    END IF

    ! some pre-calculations for the pressure terms in the u- and v-equation
    !-----------------------------------------------------------------------

    IF ( (n_step == 1) .AND. l_no_div_damping_in_1st_step ) THEN

      DO k=1, ke
        DO j=1, je
          DO i=1, ie
            ! (not beautiful, but most simple for programming:)
            pp_minus_div(i,j,k) = pp(i,j,k)  
          END DO
        END DO
      END DO

    ELSE

      DO k=1, ke
        DO j=1, je
          DO i=1, ie
            pp_minus_div(i,j,k) = pp(i,j,k)          &
                      - alpha_div_h(i,j,k) * rho(i,j,k) * div_expl_uv(i,j,k)
          END DO
        END DO
      END DO

    END IF

    IF ( .NOT. l_hor_p_grad_Mahrer ) THEN

      ! vertical derivative of the pressure perturbation:
      ! (is only needed for k=vcoord%kflat...ke)

      DO j=j_start_v, jendv + 1

        DO k=MAX(2,vcoord%kflat), ke
          DO i=i_start_u, iendu + 1
            ! work3 = p' at w-position
            work3(i,k) =              wgtfac(i,j,k)  * pp_minus_div(i,j,k)   &
                       +  (1.0_ireals-wgtfac(i,j,k)) * pp_minus_div(i,j,k-1)
          END DO
        END DO

        DO k=MAX(2,vcoord%kflat), ke-1
          DO i=i_start_u, iendu + 1
            ! dp'/dz at the scalar position (i,j,k)
            dpp_dz_s(i,j,k) = -sqrtg_r_s(i,j,k) * (work3(i,k+1) - work3(i,k) )
          END DO
        END DO

      END DO

      IF ( vcoord%kflat==1 ) THEN
        ! 1st order one-sided finite difference at top boundary (k=1):
        DO j=j_start_v, jendv + 1
          DO i=i_start_u, iendu + 1
            dpp_dz_s(i,j,1) = -sqrtg_r_s(i,j,1) * ( pp_minus_div(i,j,1)   &
                                                  - pp_minus_div(i,j,2) )
          END DO
        END DO
      ELSE
        ! dpp_dz_s(k=1) not needed, if vcoord%kflat > 1
      END IF

      IF ( .NOT. ldyn_bbc ) THEN
        ! 2nd order one-sided finite difference at bottom boundary (k=ke):
        DO j=j_start_v, jendv + 1
          DO i=i_start_u, iendu + 1
            dpp_dz_s(i,j,ke) =                                        &
                  ( wgt_one_sided(i,j,1) * pp_minus_div(i,j,ke)       &
                  + wgt_one_sided(i,j,2) * pp_minus_div(i,j,ke-1)     &
                  + wgt_one_sided(i,j,3) * pp_minus_div(i,j,ke-2) )
          END DO
        END DO
      ELSE
        ! see  'CALL dynamic_bottom_BC_for_p'  below
      END IF

    ELSE   !  l_hor_p_grad_Mahrer==.TRUE.:

      ! calculate   dp'/dz * ( z(ke-1) - z(ke) )   in k=ke:

      IF ( .NOT. ldyn_bbc ) THEN
        ! 2nd order one-sided finite difference at bottom boundary (k=ke):
        DO j=1, je
          DO i=1, ie

            dpp_dz_s(i,j,ke) =                                        &
                  ( wgt_one_sided(i,j,1) * pp_minus_div(i,j,ke)       &
                  + wgt_one_sided(i,j,2) * pp_minus_div(i,j,ke-1)     &
                  + wgt_one_sided(i,j,3) * pp_minus_div(i,j,ke-2) )

            dpp_dz_s(i,j,ke) = dpp_dz_s(i,j,ke) *    &
                  0.5_ireals * ( hhl(i,j,ke-1) - hhl(i,j,ke+1) )

          END DO
        END DO
      ELSE
        ! see  'CALL dynamic_bottom_BC_for_p'  below
      END IF

    END IF

    !-------------------------------------------------------------------------
    ! section: calculate  u^(n+1) and v^(n+1)
    !-------------------------------------------------------------------------
    ! = there is only an 'explicit' part of the u-/v-equation
    ! = forward calculation, i.e. before the implicit solution for w):

    IF ( l_no_div_damping_in_1st_step .AND. (n_step == 1) ) THEN
      switch_div_damp = 0.0_ireals
    ELSE 
      switch_div_damp = 1.0_ireals
    END IF

    IF ( ldyn_bbc ) THEN
      k_bottom = ke-1
    ELSE 
      k_bottom = ke
    END IF

    ! The calculation of the horizontal gradients of p' (and div) takes 
    ! place in 2 (or 3) layer types:
    ! (this is partly motivated by efficiency)

    ! Layer type a):  1 <= k <= k1: 
    !   no metric terms necessary, and consequently 
    !   no Mahrer discretization necessary
    k1 = MIN( vcoord%kflat-1, k_bottom )

    ! Layer type b):  k1+1 <= k <= k2
    !   metric terms are needed, but Mahrer discretization is 
    !   not necessary (it has no advantage towards the conventional discret.)
    ! Layer type c):  k2+1 <= k <= k_bottom
    !   case 'no Mahrer': same as layer type b)
    !   case 'Mahrer'   : do Mahrer-discretization

    IF ( l_hor_p_grad_Mahrer ) THEN
      IF ( l_Mahrer_opti ) THEN
        ! with this version max w is much higher than without the intermediate layer:
        k2 = MIN( k_index_same_level_Mahrer, k_bottom )
      ELSE
        ! you are on the safe side with
        k2 = k1    ! i.e. the Mahrer-discretization immediately starts below vcoord%kflat
      END IF
    ELSE
      k2 = k_bottom
    END IF


    ! --- Layer type a) ---

    DO k=1, k1

      DO j=jstart_1, jend_1
        DO i=istart_1, iend_1

          ! --- u-equation ----

          dp_dx_rho = frac_1_rhodx(i,j,k) *                    &
            &  (pp_minus_div(i+1,j,k) - pp_minus_div(i,j,k) ) ! *beta_s_1_expl(=1)

          contrib_u(i,j) = du_dt_slow(i,j,k) - dp_dx_rho

          ! --- v-equation ---

          dp_dy_rho = frac_1_rhody(i,j,k) *                    &
            &  (pp_minus_div(i,j+1,k) - pp_minus_div(i,j,k))  ! *beta_s_1_expl(=1)

          contrib_v(i,j) = dv_dt_slow(i,j,k) - dp_dy_rho

        END DO
      END DO

      DO j=jstartu, jendu
        DO i=i_start_u, iendu
          u(i,j,k) = u(i,j,k) + dts * contrib_u(i,j)
        END DO
      END DO

      DO j=j_start_v, jendv
        DO i=istartv, iendv
          v(i,j,k) = v(i,j,k) + dts * contrib_v(i,j)
        END DO
      END DO

    END DO

    ! --- Layer type  b) (and eventually c) ) ---

    DO k=k1+1, k2

      DO j=jstart_1, jend_1
        DO i=istart_1, iend_1

          ! --- u-equation ----

          ! 2 * dpp/dz  at u-position (i+1/2,j,ke):
          dpp_dz_u_double = dpp_dz_s(i,j,k) + dpp_dz_s(i+1,j,k)

          dpp_dlam_metr = dz_lam_u_half(i,j,k) * dpp_dz_u_double

          dp_dx_rho = frac_1_rhodx(i,j,k) * (                     &
                ( pp_minus_div(i+1,j,k) - pp_minus_div(i,j,k) )   & ! * beta_s_1_expl (=1)
                         - dpp_dlam_metr )                          ! * beta_s_7_expl (=1)

          contrib_u(i,j) = du_dt_slow(i,j,k) - dp_dx_rho

          ! --- v-equation ---

          ! 2 * dpp/dz  at v-position (i,j+1/2,ke):
          dpp_dz_v_double = dpp_dz_s(i,j,k) + dpp_dz_s(i,j+1,k)

          dpp_dphi_metr = dz_phi_v_half(i,j,k) * dpp_dz_v_double

          dp_dy_rho = frac_1_rhody(i,j,k) * (                     &
                ( pp_minus_div(i,j+1,k) - pp_minus_div(i,j,k) )   & ! * beta_s_1_expl (=1)
                         - dpp_dphi_metr )                          ! * beta_s_7_expl (=1)

          contrib_v(i,j) = dv_dt_slow(i,j,k) - dp_dy_rho

        END DO
      END DO

      DO j=jstartu, jendu
        DO i=i_start_u, iendu
          u(i,j,k) = u(i,j,k) + dts * contrib_u(i,j)
        END DO
      END DO

      DO j=j_start_v, jendv
        DO i=istartv, iendv
          v(i,j,k) = v(i,j,k) + dts * contrib_v(i,j)
        END DO
      END DO

    END DO

    ! --- Layer type c) ---

    IF ( l_hor_p_grad_Mahrer ) THEN

      DO k=k2+1, k_bottom

        DO j=jstart_1, jend_1
          DO i=istart_1, iend_1

            IF ( .TRUE. ) THEN
              !IF ( k <= k_star_u(i,j) ) THEN
              ! Interpolation or weak extrapolation of p 
              ! (or p-alpha rho div) on both sides:

              m = MIN( p_hor_idx(1,i,j,k), ke )
              p_r =             p_hor_wght(1,i,j,k) *pp_minus_div(i+1,j,m-1) &
                  + (1.0_ireals-p_hor_wght(1,i,j,k))*pp_minus_div(i+1,j,m)

              m = MIN( p_hor_idx(2,i,j,k), ke )
              p_l =             p_hor_wght(2,i,j,k) *pp_minus_div(i  ,j,m-1) &
                  + (1.0_ireals-p_hor_wght(2,i,j,k))*pp_minus_div(i  ,j,m)

              dp_dx_rho = frac_1_rhodx(i,j,k) * ( p_r - p_l )

            ELSE
              zu        = 0.5_ireals * ( hhl(i,j,k) + hhl(i+1,j,k) )
              zu_k_star = 0.5_ireals * ( hhl(i,j,k_star_u(i,j)) + hhl(i+1,j,k_star_u(i,j)) )
              dp_dx_rho = frac_1_rhodx(i,j,k) * ( dp_dlam_k_star(i,j)    &
                              + dbuoy_dlam_k_star(i,j) * ( zu - zu_k_star ) )

            END IF
            contrib_u(i,j) = du_dt_slow(i,j,k) - dp_dx_rho

            IF ( .TRUE. ) THEN
              !IF ( k <= k_star_v(i,j) ) THEN
              ! Interpolation or weak extrapolation of p 
              ! (or p-alpha rho div) on both sides:

              m = MIN( p_hor_idx(3,i,j,k), ke )
              p_r =             p_hor_wght(3,i,j,k) *pp_minus_div(i,j+1,m-1) &
                  + (1.0_ireals-p_hor_wght(3,i,j,k))*pp_minus_div(i,j+1,m)

              m = MIN( p_hor_idx(4,i,j,k), ke )
              p_l =             p_hor_wght(4,i,j,k) *pp_minus_div(i,j  ,m-1) &
                  + (1.0_ireals-p_hor_wght(4,i,j,k))*pp_minus_div(i,j  ,m)

              dp_dy_rho = frac_1_rhody(i,j,k) * ( p_r - p_l )

            ELSE
              zv        = 0.5_ireals * ( hhl(i,j,k) + hhl(i,j+1,k) )
              zv_k_star = 0.5_ireals * ( hhl(i,j,k_star_v(i,j)) + hhl(i,j+1,k_star_v(i,j)) )
              dp_dy_rho = frac_1_rhody(i,j,k) * ( dp_dphi_k_star(i,j)    &
                              + dbuoy_dphi_k_star(i,j) * ( zv - zv_k_star ) )

            END IF
            contrib_v(i,j) = dv_dt_slow(i,j,k) - dp_dy_rho

          END DO
        END DO

        DO j=jstartu, jendu
          DO i=i_start_u, iendu
            u(i,j,k) = u(i,j,k) + dts * contrib_u(i,j)
          END DO
        END DO

        DO j=j_start_v, jendv
          DO i=istartv, iendv
            v(i,j,k) = v(i,j,k) + dts * contrib_v(i,j)
          END DO
        END DO

      END DO

    END IF   ! if (l_hor_p_grad_Mahrer)

    IF ( ldyn_bbc ) THEN
      ! deliver values for u(:,:,ke) and v(:,:,ke):
      CALL dynamic_bottom_BC_for_p( pp_minus_div, 0.0_ireals )
    END IF


    IF ( l2dim ) THEN
      DO j=1, nboundlines
        u(:,jstartu-j,:) = u(:,jstartu,:)
        u(:,jendu  +j,:) = u(:,jstartu,:)
        v(:,jstartv-j,:) = v(:,jstartv,:)
        v(:,jendv  +j,:) = v(:,jstartv,:)
      END DO
      v(:,jendv,:) = v(:,jstartv,:)   ! necessary?
    END IF

    !-------------------------------------------------------------------------
    ! calculate the lhs of the implicit equation for w
    !-------------------------------------------------------------------------

    IF ( .NOT. l_calc_lhs_at_1st_RKstep ) THEN
      IF ( .NOT. l_small_pert_in_pT ) THEN
        IF ( l_no_div_damping_in_1st_step .AND. (n_step==1) ) THEN
          CALL lhs_of_tridiag_system_for_w( dts, .FALSE. )
        ELSE
          CALL lhs_of_tridiag_system_for_w( dts, l_3D_div_damping )
        END IF
      ELSE
        IF ( l_no_div_damping_in_1st_step ) THEN
          IF ( l_3D_div_damping ) THEN
            IF ( n_step == 1 ) THEN
              CALL lhs_of_tridiag_system_for_w( dts, .FALSE. )
            ELSE IF ( n_step == 2 ) THEN
              CALL lhs_of_tridiag_system_for_w( dts, .TRUE. )
            ELSE
              ! for n_step>2 use the coefficients calculated for n_step==2
            END IF
          ELSE
            IF ( n_step==1 ) THEN
              CALL lhs_of_tridiag_system_for_w( dts, .FALSE. )
            ELSE
              ! for n_step>1 use the coefficients calculated for n_step==1
            END IF
          END IF
        ELSE
          IF ( n_step==1 ) THEN
            CALL lhs_of_tridiag_system_for_w( dts, l_3D_div_damping )
          ELSE
            ! for n_step>1 use the coefficients calculated for n_step==1
          END IF
        END IF
      END IF
    ELSE
      ! This is the standard case: the subroutine is
      ! already called at the beginning of the whole RK-scheme.
    END IF

    !-------------------------------------------------------------------------
    ! calculate the rhs of the implicit equation for w
    !-------------------------------------------------------------------------


    ! 'horizontal' part of the divergence using u^(n+1), v^(n+1)
    !---------------------------------------------------------------

    ! 'horizontal' contribution to dzeta/dt from u and v 
    ! at the (small) timestep n+1:
    CALL calc_uv_surface( u, v, i_bbc_w_extrapol_uv_order, u_sfc, v_sfc )
    CALL calc_Z_horiz( u, v, u_sfc, v_sfc, Z_horiz)

    DO j=jstart, jend
      DO k=vcoord%kflat, ke
        DO i=istart, iend
          div_tilde_vertic_uv(i,j,k) = Z_horiz(i,j,k+1) - Z_horiz(i,j,k)
        END DO
      END DO
    END DO
    ! remark: to begin with, the approximated values 
    ! of Z_horiz(:,:,ke+1) are used.

    DO k=1, ke
      DO j=jstart, jend
        DO i=istart, iend

          ! necessary terms for the calculation of the divergence (div v)

          dug_dlam  = u(i  ,j,k) * sqrtg_eddlon_u(i,  j,k)           &
                    - u(i-1,j,k) * sqrtg_eddlon_u(i-1,j,k)

          dvg_dphi  = v(i,j  ,k) * sqrtg_cosphi_eddlat_v(i,j  ,k)    &  
                    - v(i,j-1,k) * sqrtg_cosphi_eddlat_v(i,j-1,k)

          div_tilde_horiz(i,j,k) = dug_dlam + dvg_dphi

        END DO
      END DO
    END DO

    IF ( l2dim ) THEN
      DO j=1, nboundlines
        div_tilde_horiz    (:,jstart-j,:) = div_tilde_horiz    (:,jstart,:)
        div_tilde_horiz    (:,jend  +j,:) = div_tilde_horiz    (:,jstart,:)

        div_tilde_vertic_uv(:,jstart-j,:) = div_tilde_vertic_uv(:,jstart,:)
        div_tilde_vertic_uv(:,jend  +j,:) = div_tilde_vertic_uv(:,jstart,:)
      END DO
    END IF

    IF ( n_step == 1 ) THEN
      ! vertical average of w at the OLD timestep
      ! (this is needed at the end again)
      DO k=1, ke
        DO j=jstart, jend
          DO i=istart, iend
            w_k_avg_old(i,j,k)  = 0.5_ireals * ( w(i,j,k+1) + w(i,j,k) ) 
          END DO
        END DO
      END DO
    END IF


    DO j=jstart, jend

      ! auxiliary fields:
      DO k=1, ke
        DO i=istart, iend

          ! parts of the implicit part of the divergence (of the p-/T-equation)
          div_part_p_T = sqrtgtot_r(i,j,k)                     &
                  * (                 div_tilde_horiz(i,j,k)   & ! *beta_s_3_impl (=1)
                    + beta_s_8_impl * div_tilde_vertic_uv(i,j,k) )
          ! remark: here an estimation for Z_horiz(:,:,ke+1)must be used,
          ! i.e. dzeta/dt=0 cannot be used.

          ! if  beta_s_9 /= beta_s_8, then you must distinguish between
          !    div_part_p and div_part_T:

          !div_part_p = sqrtgtot_r(i,j,k)                       &
          !        * (                 div_tilde_horiz(i,j,k)   & ! *beta_s_3_impl (=1)
          !          + beta_s_8_impl * div_tilde_vertic_uv(i,j,k) )

          !div_part_T = sqrtgtot_r(i,j,k)                       &
          !        * (                 div_tilde_horiz(i,j,k)   & ! *beta_s_5_impl (=1)
          !          + beta_s_9_impl * div_tilde_vertic_uv(i,j,k) )

          ! explicit parts of the pp-equation (at scalar position):
          p_eq_expl = pp(i,j,k) * dts_inv                                &
                 - gamma_p_tot(i,j,k) * div_expl_p_T(i,j,k)              &
                 + beta_b_3_expl * g * rho0(i,j,k) * w_k_avg_old(i,j,k)  &
                 + dpp_dt_slow(i,j,k)           

          ! explicit parts of the Tp-equation (at scalar position):
          T_eq_expl = Tp(i,j,k) * dts_inv                                &
                 - R_by_cv * T_tot(i,j,k) * div_expl_p_T(i,j,k)          &
                 - beta_b_4_expl * dT0dz(i,j,k) * w_k_avg_old(i,j,k)     &
                 + dTp_dt_slow(i,j,k)         

          ! if  beta_s_6 /= beta_s_4, then you must distinguish between
          !   div_expl_p and div_expl_T:

          ! explicit parts of the pp-equation (at scalar position):
          !p_eq_expl = pp(i,j,k) * dts_inv                                &
          !       - gamma_p_tot(i,j,k) * div_expl_p(i,j,k)                &
          !       + beta_b_3_expl * g * rho0(i,j,k) * w_k_avg_old(i,j,k)  &
          !       + dpp_dt_slow(i,j,k)           

          ! explicit parts of the Tp-equation (at scalar position):
          !T_eq_expl = Tp(i,j,k) * dts_inv                                &
          !       - R_by_cv * T_tot(i,j,k) * div_expl_T(i,j,k)            &
          !       - beta_b_4_expl * dT0dz(i,j,k) * w_k_avg_old(i,j,k)     &
          !       + dTp_dt_slow(i,j,k)         

          ! auxiliary field
          !   b_4p  :=  b_p - ( A_pu * u^(n+1) + A_pv * v^(n+1) )
          ! (is needed later for the calculation of p')

          b_4p(i,j,k) = p_eq_expl - gamma_p_tot(i,j,k) * div_part_p_T
          !b_4p(i,j,k) = p_eq_expl - gamma_p_tot(i,j,k) * div_part_p

          ! auxiliary field 
          ! work3 = [ b_T - ( A_Tu * u^(n+1) + A_Tv * v^(n+1) ) ] :

          work3(i,k) = T_eq_expl - R_by_cv * div_part_p_T * T_tot(i,j,k)
          !work3(i,k) = T_eq_expl - R_by_cv * div_part_T * T_tot(i,j,k)

        END DO
      END DO


      DO k=2, ke
        !CDIR ON_ADB(work3)
        !CDIR ON_ADB(pp)
        DO i=istart, iend

          one_m_wgtfac = 1.0_ireals - wgtfac(i,j,k)

          ! explicit parts of the w-equation (at w-position):

          ! T', p' at w-pos.:
          Tp_k_avg = wgtfac(i,j,k) * Tp(i,j,k)    &
                    + one_m_wgtfac * Tp(i,j,k-1)

          pp_k_avg = wgtfac(i,j,k) * pp(i,j,k)    &
                    + one_m_wgtfac * pp(i,j,k-1)

          ! buoyancy:
          buoy = g * ( beta_b_1_expl * p0_pT0_k_avg(i,j,k) * Tp_k_avg      &
                     - beta_b_2_expl * p_inv_k_avg(i,j,k) * pp_k_avg       &
                      + rho0_rho_q_x_k_avg(i,j,k)  )

          b_w =   w(i,j,k) * dts_inv                                       &
                + rho_inv_k_avg(i,j,k) * beta_s_2_expl * sqrtg_r_w(i,j,k)  &
                               * ( pp(i,j,k) - pp(i,j,k-1) )               &
                + buoy                                                     &
                + dw_dt_slow(i,j,k)

          ! contributions from -A_wp ( b_p - (A_pu*u^(n+1) + A_pv*v^(n+1) ) ):
          contr_b_4p =                                                     &
                 rho_inv_k_avg(i,j,k) * sqrtg_r_w(i,j,k) * beta_s_2_impl   &
                                 * ( b_4p(i,j,k) - b_4p(i,j,k-1) )         &
               - g * p_inv_k_avg(i,j,k) * beta_b_2_impl                    &
                           * ( wgtfac(i,j,k) * b_4p(i,j,k)                 &
                             +  one_m_wgtfac * b_4p(i,j,k-1) )

          ! contributions from -A_wT ( b_T - (A_Tu*u^(n+1) + A_Tv*v^(n+1) ) ):
          contr_b_T =                                             &
                 g * beta_b_1_impl * p0_pT0_k_avg(i,j,k)          &
                        * (  wgtfac(i,j,k) * work3(i,k)           &
                            + one_m_wgtfac * work3(i,k-1) )

          lgs_rhs(i,j,k) = b_w * dts_inv  + contr_b_4p + contr_b_T

        END DO
      END DO

      IF ( l_3D_div_damping .AND.                                              &
             (      (       l_no_div_damping_in_1st_step .AND. ( n_step>1 ) )  &
               .OR. ( .NOT. l_no_div_damping_in_1st_step )    )                &
         ) THEN

        ! contributions from  - 1/dts * ( A_wu u^(n+1) + A_wv v^(n+1) )
        !       and parts of  - 1/dts * b_w    :
        DO k=1, ke
          DO i=istart, iend

            work3(i,k) = alpha_div_v_to_h * alpha_div_h(i,j,k) * rho(i,j,k)  &
                  * ( sqrtgtot_r(i,j,k)                                      &
                        * ( beta_d_3_impl * div_tilde_horiz(i,j,k)           &
                          + beta_d_6_impl * div_tilde_vertic_uv(i,j,k) )     &
                    + div_expl_w(i,j,k)    )
          END DO
        END DO

        DO k=2, ke
          DO i=istart, iend

            lgs_rhs(i,j,k) = lgs_rhs(i,j,k)                              &
                   - dts_inv * rho_inv_k_avg(i,j,k) * sqrtg_r_w(i,j,k)   &
                         * ( work3(i,k) - work3(i,k-1) )

          END DO
        END DO

      END IF

    END DO

    ! --- at top boundary w=0: ---
    lgs_rhs(:,:,1) = 0.0_ireals

    ! --- at bottom boundary (d zeta/dt=0): ---
    CALL bottom_BC_for_w( u_sfc, v_sfc, dts, lgs_rhs(:,:,ke1) )


    ! solve the implicit equation for w
    !---------------------------------------------------------------

    IF ( l_use_tridag ) THEN
      CALL tridag( lgs_A, lgs_B, lgs_C, lgs_rhs, w, ke1 )
    ELSE

      ! 2nd part of the LU decomposition: solve the bidiagonal equation for w
      ! lgs_rhs(i,j,ke+1)=w(i,j,ke+1) is set in Subr. 'bottom_BC_for_w'

      DO k= ke, 2, -1
        DO j=jstart, jend
          DO i=istart, iend
            lgs_rhs(i,j,k) = (lgs_rhs(i,j,k) - lgs_C(i,j,k) * lgs_rhs(i,j,k+1) ) &
                      * lgs_B(i,j,k)     
          END DO
        END DO
      END DO

      ! upper boundary condition:
      ! w(:,:,1)    = 0.0_ireals
      w( istart:iend, jstart:jend, 1) = 0.0_ireals

      DO k= 2, ke
        DO j=jstart, jend
          DO i=istart, iend
            w(i,j,k) = lgs_rhs(i,j,k) - lgs_A(i,j,k) * w(i,j,k-1)
          END DO
        END DO
      END DO
      ! lower boundary condition:
      ! w(:,:,ke+1) = lgs_rhs(:,:,ke+1)
      w(istart:iend,jstart:jend, ke+1) =lgs_rhs(istart:iend,jstart:jend, ke+1)

    END IF

    ! Rayleigh damping on w according to Klemp et al. (2008; MWR)
    IF (lspubc .AND. (itype_spubc == 3) ) THEN
      DO k = 1, ke_rd
        DO  j = jstart, jend
          DO  i = istart, iend
            ! wrdfac is defined on top for efficiency reasons
            w(i,j,k) = w(i,j,k)*wrdfac(k)  
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF ( l2dim ) THEN
      DO j=1, nboundlines
        w(:,jstart-j,:) = w(:,jstart,:)
        w(:,jend  +j,:) = w(:,jstart,:)
      END DO
    END IF


    ! calculate new fields pp and Tp
    !---------------------------------------------------------------

    IF ( i_div_at_horiz_bound==1 ) THEN
      k_bottom = ke-1
    ELSE
      k_bottom = ke
    END IF

    DO k=1, k_bottom
      DO j=jstart, jend
        DO i=istart, iend

          ! vertical average of w at the new time step:
          w_k_avg_new = 0.5_ireals * ( w(i,j,k+1) + w(i,j,k) )

          ! --- p-equation --------------

          dw_dzeta = w(i,j,k+1) - w(i,j,k)

          pp(i,j,k) = dts * ( b_4p(i,j,k)                               &
                    + gamma_p_tot(i,j,k) * beta_s_4_impl                &
                            * sqrtg_r_s(i,j,k) * dw_dzeta               &
                    + beta_b_3_impl * g * rho0(i,j,k) * w_k_avg_new )

          ! --- T-equation --------------

          div_tilde_vertic_w  = -r_cosphi(j) * dw_dzeta

          div_impl_T = sqrtgtot_r(i,j,k) * (                 &
                                 div_tilde_horiz(i,j,k)      & ! * beta_s_5_impl (=1)
               + beta_s_9_impl * div_tilde_vertic_uv(i,j,k)  &
               + beta_s_6_impl * div_tilde_vertic_w ) 

          Tp(i,j,k) = Tp(i,j,k) + dts *                                      &
              (-R_by_cv * T_tot(i,j,k) * (div_expl_p_T(i,j,k) + div_impl_T ) &
               -dT0dz(i,j,k) * ( beta_b_4_impl * w_k_avg_new                 &
                                 + beta_b_4_expl * w_k_avg_old(i,j,k) )      &
                + dTp_dt_slow(i,j,k) )

          ! in the case beta_s_6 /= beta_s_4 distinguish between 
          !   div_expl_p and div_expl_T:
          !Tp(i,j,k) = Tp(i,j,k) + dts *                                     &
          !    (-R_by_cv * T_tot(i,j,k) * ( div_expl_T(i,j,k) + div_impl_T ) &
          !     -dT0dz(i,j,k) * ( beta_b_4_impl * w_k_avg_new                &
          !                       + beta_b_4_expl * w_k_avg_old(i,j,k) )     &
          !      + dTp_dt_slow(i,j,k) )

          ! --- save avg(w) for the next time step: ---
          w_k_avg_old(i,j,k) = w_k_avg_new

        END DO
      END DO
    END DO

    IF ( i_div_at_horiz_bound==1 ) THEN
      ! extra treatment at the lower boundary:
      k = ke

      DO j=jstart, jend
        DO i=istart, iend

          ! vertical average of w at the new time step:
          w_k_avg_new = 0.5_ireals * ( w(i,j,k+1) + w(i,j,k) )

          ! --- p-equation --------------

          !dw_dzeta = w(i,j,k+1) - w(i,j,k)

          ! version 1 (use of exactly dzeta/dt=0):
          div_tilde_vertic  = - ( Z_horiz(i,j,ke) - r_cosphi(j) * w(i,j,ke) )
          div_tilde_vertic_w = div_tilde_vertic - div_tilde_vertic_uv(i,j,ke)

          pp(i,j,k) = dts * ( b_4p(i,j,k)                               &
                    - gamma_p_tot(i,j,k) * beta_s_4_impl                &
                        * sqrtgtot_r(i,j,ke) * div_tilde_vertic_w       &
                    + beta_b_3_impl * g * rho0(i,j,k) * w_k_avg_new )

          ! --- T-equation --------------

          div_impl_T = sqrtgtot_r(i,j,k) * (            &
                                 div_tilde_horiz(i,j,k) & ! *beta_s_5_impl(=1)
               + beta_s_6_impl * div_tilde_vertic ) 

          Tp(i,j,k) = Tp(i,j,k) + dts *                                     &
              (-R_by_cv * T_tot(i,j,k) * (div_expl_p_T(i,j,k) + div_impl_T) &
               -dT0dz(i,j,k) * ( beta_b_4_impl * w_k_avg_new                &
                                 + beta_b_4_expl * w_k_avg_old(i,j,k) )     &
                + dTp_dt_slow(i,j,k) )

          ! in the case beta_s_6 /= beta_s_4 distinguish between 
          !   div_expl_p and div_expl_T:
          !Tp(i,j,k) = Tp(i,j,k) + dts *                                    &
          !    (-R_by_cv * T_tot(i,j,k) * (div_expl_T(i,j,k) + div_impl_T)  &
          !     -dT0dz(i,j,k) * ( beta_b_4_impl * w_k_avg_new               &
          !                       + beta_b_4_expl * w_k_avg_old(i,j,k) )    &
          !      + dTp_dt_slow(i,j,k) )

          ! --- save avg(w) for the next time step: ---
          w_k_avg_old(i,j,k) = w_k_avg_new

        END DO
      END DO
    END IF

    IF ( l2dim ) THEN
      DO j=1, nboundlines
        pp(:,jstart-j,:) = pp(:,jstart,:)
        pp(:,jend  +j,:) = pp(:,jstart,:)
        Tp(:,jstart-j,:) = Tp(:,jstart,:)
        Tp(:,jend  +j,:) = Tp(:,jstart,:)
      END DO
    END IF

    IF ( lradlbc .AND. .NOT.( lperi_x .OR. lperi_y ) ) THEN
      ! radiative lateral boundary condition (init)
      ! update of the boundary fields
      CALL set_radiative_lbc( u, v, w, pp, Tp,                         &
                              bd_west, bd_east, bd_north, bd_south )
    END IF

  END DO   ! loop: n_step

  !---------------------------------------------------------------
  ! end of the time integration
  !---------------------------------------------------------------


  !---------------------------------------------------------------
  ! Set free-slip lateral boundary conditions on w for
  ! non-periodic boundary conditions
  ! (applies also for single-processor and/or 2D case)
  !---------------------------------------------------------------

  IF (lw_freeslip .AND. .NOT. lperi_x ) THEN

    IF (my_cart_neigh(1) == -1) THEN
      DO  k = 1, ke1
        DO i = 1, nboundlines
!CDIR NOLOOPCHG
          DO  j = jstart, jend
            w(i,j,k) = w(istart,j,k)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (my_cart_neigh(3) == -1) THEN
      DO  k = 1, ke1
        DO i = ie-nboundlines+1, ie
!CDIR NOLOOPCHG
          DO  j = jstart, jend
            w(i,j,k) = w(iend,j,k)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  ENDIF

  IF ( lw_freeslip .AND. .NOT. lperi_y ) THEN

    IF (my_cart_neigh(4) == -1) THEN
      DO  k = 1, ke1
        DO  j = 1, nboundlines
          w(:,j,k) = w(:,jstart,k)
        ENDDO
      ENDDO
    ENDIF

    IF (my_cart_neigh(2) == -1) THEN
      DO  k = 1, ke1
        DO  j = je-nboundlines+1, je
          w(:,j,k) = w(:,jend,k)
        ENDDO
      ENDDO
    ENDIF

  ENDIF


  ! deallocate fields
  !---------------------------------------------------------------

  DEALLOCATE( lgs_rhs )

  IF ( .NOT. l_calc_lhs_at_1st_RKstep ) THEN
    CALL dealloc_fast_waves_sc 
  END IF

  DEALLOCATE ( u_sfc )
  DEALLOCATE ( v_sfc )

  IF ( l_3D_div_damping ) THEN
    DEALLOCATE( div_expl_w   )
  END IF

  IF ( l_hor_p_grad_Mahrer ) THEN
    !DEALLOCATE( hml )
    !DEALLOCATE( hml_12_inv )
    !DEALLOCATE( hml_13_inv )
    !DEALLOCATE( hml_23_inv )
  END IF

  DEALLOCATE( dpp_dz_s )

  IF ( l_hor_p_grad_Mahrer ) THEN
    DEALLOCATE( dp_dlam_k_star   , STAT=izstata)
    DEALLOCATE( dp_dphi_k_star   , STAT=izstata)
    DEALLOCATE( dbuoy_dlam_k_star, STAT=izstata)
    DEALLOCATE( dbuoy_dphi_k_star, STAT=izstata)
  END IF


  IF ( lradlbc ) THEN
    DEALLOCATE (  bd_west )
    DEALLOCATE (  bd_east )
    DEALLOCATE (  bd_north )
    DEALLOCATE (  bd_south )
  END IF

  !--------------------------------------------------------------------------
  !  End of the module procedure fast_waves_strong_conserv
  !--------------------------------------------------------------------------

CONTAINS

  !==========================================================================

  SUBROUTINE dynamic_bottom_BC_for_p( pp_mod, switch_div_damp )

    !------------------------------------------------------------------------
    ! Description:
    !
    !   boundary treatment for d p'/dz at the lower boundary
    !
    ! Input:
    !   pp(:,:,:)
    !   switch_div_damp
    !   others:
    !
    ! Output:
    !   u(:,:,ke), v(:,:,ke)
    !
    ! Methods:
    !   see A. Gassmann, COSMO-Newsletter
    !
    !------------------------------------------------------------------------

    IMPLICIT NONE

    REAL(KIND=ireals), INTENT(in) :: pp_mod(ie,je,ke)
    REAL(KIND=ireals), INTENT(in) :: switch_div_damp 

    REAL (KIND=ireals), ALLOCATABLE ::               &
      &   denom_u(:,:), denom_v(:,:),                &
      &   contr_u(:,:), contr_v(:,:), contr_w(:,:)

    INTEGER :: i, j 
    REAL (KIND=ireals) :: dz_dlam_u,  dz_dphi_u
    REAL (KIND=ireals) :: dz_dlam_v,  dz_dphi_v
    REAL (KIND=ireals) :: rho_inv_u, rho_inv_v
    REAL (KIND=ireals) :: dp_dx, dp_dy
    REAL (KIND=ireals) :: au, av
    REAL (KIND=ireals) :: frac_1_8
    REAL (KIND=ireals) :: dpp_dz_u, dpp_dz_v

    ALLOCATE( denom_u (ie, je) )
    ALLOCATE( denom_v (ie, je) )
    ALLOCATE( contr_u (ie, je) )
    ALLOCATE( contr_v (ie, je) )
    ALLOCATE( contr_w (ie, je) )

    frac_1_8 = 1.0_ireals / 8.0_ireals


    DO j=jstartu, jendu
      DO i=i_start_u, iendu

        dz_dlam_u = 0.5_ireals * ( dz_dlam(i,j,ke)   + dz_dlam(i,j,ke1) )

        dz_dphi_u = ( dz_dphi(i,j  ,ke)  + dz_dphi(i+1,j  ,ke)       &
                    + dz_dphi(i,j-1,ke)  + dz_dphi(i+1,j-1,ke)       & 
                    + dz_dphi(i,j  ,ke1) + dz_dphi(i+1,j  ,ke1)      &
                    + dz_dphi(i,j-1,ke1) + dz_dphi(i+1,j-1,ke1) ) * frac_1_8

        denom_u(i,j) = ( acrlat(j,1)   * dz_dlam_u )**2   &
                     + ( r_earth_recip * dz_dphi_u )**2   &
                     + 1.0_ireals

        rho_inv_u = 0.5_ireals * ( rho_inv(i+1,j,ke) + rho_inv(i,j,ke) )

        denom_u(i,j) = denom_u(i,j) * rho_inv_u

      END DO
    END DO


    DO j=j_start_v, jendv
      DO i=istartv, iendv

        dz_dlam_v = ( dz_dlam(i,j  ,ke)  + dz_dlam(i-1,j  ,ke)     &
                    + dz_dlam(i,j+1,ke)  + dz_dlam(i-1,j+1,ke)     &
                    + dz_dlam(i,j  ,ke1) + dz_dlam(i-1,j  ,ke1)    &
                    + dz_dlam(i,j+1,ke1) + dz_dlam(i-1,j+1,ke1) ) * frac_1_8

        dz_dphi_v = 0.5_ireals * ( dz_dphi(i,j,ke)   + dz_dphi(i,j,ke1) )

        denom_v(i,j) = ( acrlat(j,2)   * dz_dlam_v )**2   &
                     + ( r_earth_recip * dz_dphi_v )**2   &
                     + 1.0_ireals

        rho_inv_v = 0.5_ireals * ( rho_inv(i,j+1,ke) + rho_inv(i,j,ke) )

        denom_v(i,j) = denom_v(i,j) * rho_inv_v

      END DO
    END DO

    DO j = MIN(jstartu, j_start_v), MAX( jendu, jendv+1)
      DO i = MIN(i_start_u, istartv-1),  MAX(iendu, iendv)
        ! contributions from the u-eqn. (at the u-position)

        rho_inv_u = 0.5_ireals * ( rho_inv(i+1,j,ke) + rho_inv(i,j,ke) )

        dp_dx = acrlat(j,1) * ( pp_mod(i+1,j,ke) - pp_mod(i,j,ke) ) * eddlon  ! * beta_s_1_expl (=1)

        contr_u(i,j) = -rho_inv_u * dp_dx + du_dt_slow(i,j,ke)

      END DO
    END DO

    DO j = MIN(j_start_v, jstartu-1), MAX(jendv, jendu)
      DO i = MIN(istartv, i_start_u), MAX(iendv, iendu+1)
        ! contributions from the v-eqn. (at the v-position)

        rho_inv_v = 0.5_ireals * ( rho_inv(i,j+1,ke) + rho_inv(i,j,ke) )

        dp_dy = r_earth_recip * ( pp_mod(i,j+1,ke) - pp_mod(i,j,ke) ) * eddlat  ! * beta_s_1_expl (=1)

        contr_v(i,j) = -rho_inv_v * dp_dy + dv_dt_slow(i,j,ke)

      END DO
    END DO

    DO j = MIN(jstartu, j_start_v), MAX( jendu, jendv+1)
      DO i = MIN(i_start_u, istartv), MAX(iendu+1, iendv)
        ! contributions from the w-eqn. (at the s-position!)

        ! remark: you cannot use Tp because it is not exchanged; therefore
        ! approximate it by T_tot from the beginning of the time step
        contr_w(i,j) = g * ( -1.0_ireals                                        &
              + p0(i,j,ke) / ( p0(i,j,ke) + pp(i,j,ke) )                        &
                 * T_tot(i,j,ke) / T0(i,j,ke)                                   &
                 * ( 1.0_ireals + ( rvd_m_o * qv(i,j,ke) - q_cond(i,j,ke) ) ) ) &
              + 0.5_ireals * ( dw_dt_slow(i,j,ke1) + dw_dt_slow(i,j,ke) )

        !   - alpha_div_v_to_h * alpha_div_h(i,j,ke) * sqrtg_r_s(i,j,ke) *   &
        !       ( div_expl_w(i,j,ke) - 0.0_ireals )             &

      END DO
    END DO

    DO j=jstartu, jendu
      DO i=i_start_u, iendu

        dz_dlam_u = 0.5_ireals * ( dz_dlam(i,j,ke1) + dz_dlam(i,j,ke) )
        au = acrlat(j,1) * dz_dlam_u * contr_u(i,j)

        av = frac_1_8 * r_earth_recip * (                                    &
               (dz_dphi(i  ,j  ,ke1)+dz_dphi(i  ,j  ,ke)) *contr_v(i  ,j)    &
             + (dz_dphi(i+1,j  ,ke1)+dz_dphi(i+1,j  ,ke)) *contr_v(i+1,j)    &
             + (dz_dphi(i  ,j-1,ke1)+dz_dphi(i  ,j-1,ke)) *contr_v(i  ,j-1)  &
             + (dz_dphi(i+1,j-1,ke1)+dz_dphi(i+1,j-1,ke)) *contr_v(i+1,j-1) )

        dpp_dz_u =  ( -au - av                                            &
                      + 0.5_ireals  * ( contr_w(i,j) + contr_w(i+1,j) ) ) &
                / denom_u(i,j)

        rho_inv_u = 0.5_ireals * ( rho_inv(i+1,j,ke) + rho_inv(i,j,ke) )

        u(i,j,ke) = u(i,j,ke) + dts * ( contr_u(i,j)  &
                          + rho_inv_u * acrlat(j,1) * dz_dlam_u * dpp_dz_u )

      END DO
    END DO

    DO j=j_start_v, jendv
      DO i=istartv, iendv

        au = frac_1_8 * (                                                                       &
              acrlat(j  ,1) *((dz_dlam(i  ,j  ,ke1) + dz_dlam(i  ,j  ,ke)) *contr_u(i  ,j)      &
                             +(dz_dlam(i-1,j  ,ke1) + dz_dlam(i-1,j  ,ke)) *contr_u(i-1,j) )    &
            + acrlat(j+1,1) *((dz_dlam(i  ,j+1,ke1) + dz_dlam(i  ,j+1,ke)) *contr_u(i  ,j+1)    &
                             +(dz_dlam(i-1,j+1,ke1) + dz_dlam(i-1,j+1,ke)) *contr_u(i-1,j+1) ) )

        dz_dphi_v = 0.5_ireals * ( dz_dphi(i,j,ke1) + dz_dphi(i,j,ke) )
        av = r_earth_recip * dz_dphi_v * contr_v(i,j)

        dpp_dz_v =  ( -au - av                                            &
                      + 0.5_ireals  * ( contr_w(i,j) + contr_w(i,j+1) ) ) &
                / denom_v(i,j)

        rho_inv_v = 0.5_ireals * ( rho_inv(i,j+1,ke) + rho_inv(i,j,ke) )

        v(i,j,ke) = v(i,j,ke) + dts * ( contr_v(i,j)         &
                         + rho_inv_v * r_earth_recip * dz_dphi_v * dpp_dz_v )

      END DO
    END DO



    DEALLOCATE( denom_u )
    DEALLOCATE( denom_v )
    DEALLOCATE( contr_u )
    DEALLOCATE( contr_v )
    DEALLOCATE( contr_w )

  END SUBROUTINE dynamic_bottom_BC_for_p

  !===========================================================================


END SUBROUTINE fast_waves_strong_conserv


!=============================================================================
!=============================================================================

SUBROUTINE calc_dz_dx( dz_dlam, dz_dphi,               &
                       dz_lam_u_half, dz_phi_v_half,   &
                       dz_dlam_4ord, dz_dphi_4ord    )

  !--------------------------------------------------------------------------
  !
  ! Description:
  !   calculate other metric terms of the terrain-following coordinate
  !   especially needed in fast_waves_strong_conserv.
  !
  !   Optionally a 4th order approximation of dz/dlam, dz_dphi
  !   is calculated additionally
  !   (only used for the free slip bottom BC for w).
  !
  ! (This subroutine would rather fit to grid_metrics_utilities.f90)
  !
  !--------------------------------------------------------------------------

  REAL (KIND=ireals), INTENT(out) :: dz_dlam       (ie,je,ke1)
  REAL (KIND=ireals), INTENT(out) :: dz_dphi       (ie,je,ke1)
  REAL (KIND=ireals), INTENT(out) :: dz_lam_u_half (ie,je,ke)
  REAL (KIND=ireals), INTENT(out) :: dz_phi_v_half (ie,je,ke)

  REAL (KIND=ireals), INTENT(out), OPTIONAL :: dz_dlam_4ord (ie,je)
  REAL (KIND=ireals), INTENT(out), OPTIONAL :: dz_dphi_4ord (ie,je)

  REAL (KIND=ireals) :: frac_1_12, frac_8_12


  REAL (KIND=ireals) :: dlambda, dphi

  INTEGER (KIND=iintegers) ::  &
    kzdims(24),          & ! Vertical dimensions for exchg_boundaries
    i,  j,  k,           & ! Loop indices
    k_top, k_bottom,     &
    izerror,             & ! error status
    izstata,             & ! error status at allocation
    izstatd                ! error status at deallocation

  CHARACTER(LEN=100) :: yzerrmsg

  IF (izdebug >= 10) THEN
    WRITE(*,*) "Subr. [calc_dz_dx] ..."
  END IF

  dz_dlam(:,:, 1:vcoord%kflat) = 0.0_ireals
  DO k=vcoord%kflat+1, ke1
    DO j=1, je
      DO i=1, ie-1

        ! J_lambda := dz / d lambda   (at uw-position (i+1/2, j, k-1/2) )
        dz_dlam(i,j,k) = ( hhl(i+1,j,k) - hhl(i,j,k) ) * eddlon

      END DO
    END DO
  END DO
  dz_dlam(ie,:,:) = 0.0_ireals   ! arbitrary value; avoids access to NaN

  dz_dphi(:,:, 1:vcoord%kflat) = 0.0_ireals
  DO k=vcoord%kflat+1, ke1
    DO j=1, je-1
      DO i=1, ie

        ! J_phi := dz / d phi         (at vw-position (i, j+1/2, k-1/2) )
        dz_dphi(i,j,k) = ( hhl(i,j+1,k) - hhl(i,j,k) ) * eddlat

      END DO
    END DO
  END DO
  dz_dphi(:,je,:) = 0.0_ireals   ! arbitrary value; avoids access to NaN


  dlambda = dlon * pi / 180.0_ireals
  dphi    = dlat * pi / 180.0_ireals

  dz_lam_u_half(:,:,1:vcoord%kflat-1) = 0.0_ireals
  dz_phi_v_half(:,:,1:vcoord%kflat-1) = 0.0_ireals
  DO k=vcoord%kflat, ke
    DO j=1, je
      DO i=1, ie
        ! dz = 1/2 * d lambda * dz/dlam  (at u position)
        dz_lam_u_half(i,j,k) = 0.25_ireals * dlambda         &
                     * ( dz_dlam(i,j,k+1) + dz_dlam(i,j,k) )

        ! dz = 1/2 * d phi    * dz/dphi  (at v position)
        dz_phi_v_half(i,j,k) = 0.25_ireals * dphi            &
                     * ( dz_dphi(i,j,k+1) + dz_dphi(i,j,k) )
      END DO
    END DO
  END DO


  IF ( l2dim ) THEN
    DO j=1, nboundlines
      dz_dlam(:,jstartu-j,:) = dz_dlam(:,jstartu,:)
      dz_dlam(:,jendu  +j,:) = dz_dlam(:,jstartu,:)

      dz_dphi(:,jstartv-j,:) = dz_dphi(:,jstartv,:)
      dz_dphi(:,jendv  +j,:) = dz_dphi(:,jstartv,:)

      dz_lam_u_half(:,jstartu-j,:) = dz_lam_u_half(:,jstartu,:)
      dz_lam_u_half(:,jendu  +j,:) = dz_lam_u_half(:,jstartu,:)

      dz_phi_v_half(:,jstartv-j,:) = dz_phi_v_half(:,jstartv,:)
      dz_phi_v_half(:,jendv  +j,:) = dz_phi_v_half(:,jstartv,:)

    END DO
    dz_dphi      (:,jendv,:) = dz_dphi      (:,jstartv,:)   ! necessary?
    dz_phi_v_half(:,jendv,:) = dz_phi_v_half(:,jstartv,:)   ! necessary?
  END IF


  IF ( PRESENT( dz_dlam_4ord ) .AND. PRESENT( dz_dphi_4ord ) ) THEN

    frac_8_12 = 8.0_ireals / 12.0_ireals
    frac_1_12 = 1.0_ireals / 12.0_ireals

    DO j=1, je
      DO i=3, ie-2
        ! J_lambda := dz / d lambda   (at w-position (i, j, ke+1/2) )
        dz_dlam_4ord(i,j) = eddlon *                              &
             ( frac_8_12 * ( hhl(i+1,j,ke1) - hhl(i-1,j,ke1) )    &
             - frac_1_12 * ( hhl(i+2,j,ke1) - hhl(i-2,j,ke1) ) )
      END DO
    END DO
    dz_dlam_4ord(1:2,    :)=0.0_ireals ! arbitrary value; avoid access to NaN
    dz_dlam_4ord(ie-1:ie,:)=0.0_ireals ! arbitrary value; avoid access to NaN

    DO j=3, je-2
      DO i=1, ie
        ! J_phi := dz / d phi         (at w-position (i, j, ke+1/2) )
        dz_dphi_4ord(i,j) = eddlat *                              &
             ( frac_8_12 * ( hhl(i,j+1,ke1) - hhl(i,j-1,ke1) )    &
             - frac_1_12 * ( hhl(i,j+2,ke1) - hhl(i,j-2,ke1) ) )
      END DO
    END DO
    dz_dphi_4ord(:,1:2)    =0.0_ireals ! arbitrary value; avoid access to NaN
    dz_dphi_4ord(:,je-1:je)=0.0_ireals ! arbitrary value; avoid access to NaN


    IF ( l2dim ) THEN
      DO j=1, nboundlines
        dz_dlam_4ord(:,jstartu-j) = dz_dlam_4ord(:,jstartu)
        dz_dlam_4ord(:,jendu  +j) = dz_dlam_4ord(:,jstartu)
        dz_dphi_4ord(:,jstartv-j) = dz_dphi_4ord(:,jstartv)
        dz_dphi_4ord(:,jendv  +j) = dz_dphi_4ord(:,jstartv)
      END DO
      dz_dphi_4ord(:,jendv) = dz_dphi_4ord(:,jstartv)   ! necessary?
    END IF

  END IF


  IF ( lperi_x .OR. lperi_y ) THEN

    IF ( PRESENT( dz_dlam_4ord ) .AND. PRESENT( dz_dphi_4ord ) ) THEN

      kzdims(1:24) =                        &
           (/  ke1, ke1 , ke, ke,   1,      &
                1,   0,   0,   0,   0,      &
                0,   0,   0,   0,   0,      &
                0,   0,   0,   0,   0,      &
                0,   0,   0,   0  /)

      ! non-optimized version (icase=0 <--> ldatatypes=.FALSE.)
      CALL exchg_boundaries                                                   &
          (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,ie,je, &
          kzdims, jstartpar, jendpar, nboundlines,nboundlines,my_cart_neigh,  &
          lperi_x, lperi_y, l2dim,                                            &
          15000+nexch_tag, .FALSE., ncomm_type, izerror, yzerrmsg,            &
          dz_dlam       (:,:,:),       &
          dz_dphi       (:,:,:),       &
          dz_lam_u_half (:,:,:),       &
          dz_phi_v_half (:,:,:),       &
          dz_dlam_4ord  (:,:),         &
          dz_dphi_4ord  (:,:)       )

    ELSE

      kzdims(1:24) =                        &
           (/  ke1, ke1 , ke, ke,   0,      &
                0,   0,   0,   0,   0,      &
                0,   0,   0,   0,   0,      &
                0,   0,   0,   0,   0,      &
                0,   0,   0,   0  /)

      ! non-optimized version (icase=0 <--> ldatatypes=.FALSE.)
      CALL exchg_boundaries                                                   &
          (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,ie,je, &
          kzdims, jstartpar, jendpar, nboundlines,nboundlines,my_cart_neigh,  &
          lperi_x, lperi_y, l2dim,                                            &
          20000+nexch_tag, .FALSE., ncomm_type, izerror, yzerrmsg,            &
          dz_dlam       (:,:,:),       &
          dz_dphi       (:,:,:),       &
          dz_lam_u_half (:,:,:),       &
          dz_phi_v_half (:,:,:) )

    END IF

    IF ( izerror /= 0_iintegers ) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'calc_dz_dx')
    END IF

  END IF

END SUBROUTINE calc_dz_dx

!=============================================================================
!=============================================================================

SUBROUTINE calc_total_T_p_rho( Tp, pp, qv, q_cond )

  !------------------------------------------------------------------
  ! Description:
  !   compute the total values of T, p, rho
  !   and some of their vertical averages
  !
  ! Input:
  !   Tp(:,:,:), pp(:,:,:), qv(:,:,:), q_cond(:,:,:)
  !
  ! Output:
  !   3D-fields (they are global in module fast_waves_strong_conserv):
  !        T_tot, gamma_p_tot, rho, rho_inv, rho_inv_k_avg,
  !        frac_1_rhodx, frac_1_rhody, 
  !        p0T_pT0_k_avg, p0_pT0_k_avg, p_inv_k_avg
  !
  !------------------------------------------------------------------

  USE data_fields, ONLY: p0hl, T0hl

  REAL (KIND=ireals), INTENT(in) :: Tp    (ie,je,ke)
  REAL (KIND=ireals), INTENT(in) :: pp    (ie,je,ke)
  REAL (KIND=ireals), INTENT(in) :: qv    (ie,je,ke)
  REAL (KIND=ireals), INTENT(in) :: q_cond(ie,je,ke)

  INTEGER :: i, j, k
  REAL (KIND=ireals) :: one_m_wgtfac
  REAL (KIND=ireals) :: r_earth_recip
  REAL (KIND=ireals) :: Tp_k_avg   ! T'   on half levels
  REAL (KIND=ireals) :: pp_k_avg   ! p'   on half levels
  REAL (KIND=ireals) :: p0_p_k_avg ! p0/p on half levels
  REAL (KIND=ireals) :: T_T0_k_avg ! T/T0 on half levels
  REAL (KIND=ireals) :: T_hl       ! T    on half levels
  REAL (KIND=ireals) :: p_hl       ! p    on half levels
  REAL (KIND=ireals) :: qv_hl      ! qv   on half levels
  REAL (KIND=ireals) :: q_cond_hl  ! condensate on half levels

  r_earth_recip = 1.0_ireals / r_earth

  IF ( izdebug >= 10 ) THEN
    WRITE(*,*) "Subr. [calc_total_T_p_rho] ..."
  END IF

  ! total temperature, pressure, and density:
  DO k=1, ke
    DO j=1,je
      DO i=1,ie
        T_tot      (i,j,k) =   T0(i,j,k) + Tp(i,j,k)
        gamma_p_tot(i,j,k) = ( p0(i,j,k) + pp(i,j,k) ) * gamma
      END DO
    END DO
  END DO

  IF ( l2dim ) THEN
    DO j=1, nboundlines
      T_tot(:,jstart-j,:) = T_tot(:,jstart,:)
      T_tot(:,jend  +j,:) = T_tot(:,jstart,:)
      gamma_p_tot(:,jstart-j,:) = gamma_p_tot(:,jstart,:)
      gamma_p_tot(:,jend  +j,:) = gamma_p_tot(:,jstart,:)
    END DO
  END IF

  CALL calrho_fw ( T_tot, gamma_p_tot, qv, q_cond, rho )

  DO k=1, ke
    DO j=1, je
      DO i=1, ie
        rho_inv(i,j,k) = 1.0_ireals / rho(i,j,k)
      END DO
    END DO
  END DO


  DO k=1, ke
    DO j=1, je-1
      DO i=1, ie-1
        ! 1 / ( rho * dx ) =
        frac_1_rhodx(i,j,k) = acrlat(j,1) * eddlon                  &
              * 0.5_ireals * ( rho_inv(i+1,j,k) + rho_inv(i,j,k) )

        ! 1 / ( rho * dy ) =
        frac_1_rhody(i,j,k) = r_earth_recip * eddlat                &
              * 0.5_ireals * ( rho_inv(i,j+1,k) + rho_inv(i,j,k) )
      END DO
    END DO
  END DO


  DO k=2, ke
    DO j=jstart, jend
      DO i=istart, iend

        one_m_wgtfac = 1.0_ireals - wgtfac(i,j,k)

        Tp_k_avg = wgtfac(i,j,k) * Tp(i,j,k)             &
                  + one_m_wgtfac * Tp(i,j,k-1)

        pp_k_avg = wgtfac(i,j,k) * pp(i,j,k)             &
                  + one_m_wgtfac * pp(i,j,k-1)

        T_hl = T0hl(i,j,k) + Tp_k_avg
        p_hl = p0hl(i,j,k) + pp_k_avg

        qv_hl = wgtfac(i,j,k) * qv(i,j,k)                &
               + one_m_wgtfac * qv(i,j,k-1)

        q_cond_hl = wgtfac(i,j,k) * q_cond(i,j,k)        &
                   + one_m_wgtfac * q_cond(i,j,k-1)

        rho_inv_k_avg(i,j,k) = 1.0_ireals /rho_fw_gp( T_hl, p_hl, qv_hl, q_cond_hl )

        p0_p_k_avg = 1.0_ireals / ( 1.0_ireals + pp_k_avg / p0hl(i,j,k) )
        T_T0_k_avg = 1.0_ireals + Tp_k_avg / T0hl(i,j,k)

        p_inv_k_avg  (i,j,k) = 1.0_ireals / p_hl
        p0_pT0_k_avg (i,j,k) = p0_p_k_avg / T0hl(i,j,k)
        p0T_pT0_k_avg(i,j,k) = p0_p_k_avg * T_T0_k_avg 

      END DO
    END DO
  END DO

END SUBROUTINE calc_total_T_p_rho

!=============================================================================
!=============================================================================

SUBROUTINE lhs_of_tridiag_system_for_w( dts, zl_3D_div_damping )

  !------------------------------------------------------------
  ! Description:
  !  determine the left hand side of the tridiagonal system
  !    A w(k-1) + B w(k) + C w(k+1) = rhs(k)
  !
  ! output:
  !   coefficient fields: lgs_A, lgs_B, lgs_C 
  !
  ! Remark: actually the coefficients a_.. must be dependent only on k.
  ! But this results in quite short vectors and additionally generates 
  ! a lot of bank conflicts on the NEC SX9. Therefore an additional 
  ! index i was introduced.
  !
  ! In the case this routine is called several times
  ! (i.e. if  l_small_pert_in_pT = .TRUE.),
  ! it is more efficient to do the allocations in 'fast_waves_strong_conserv'
  !
  !------------------------------------------------------------

  REAL (KIND=ireals), INTENT(in)  :: dts
  LOGICAL, INTENT(IN)             :: zl_3D_div_damping

  REAL (KIND=ireals), ALLOCATABLE :: a_w_p1 (:,:)
  REAL (KIND=ireals), ALLOCATABLE :: a_w_p2 (:,:)
  REAL (KIND=ireals), ALLOCATABLE :: a_p_w1 (:,:)
  REAL (KIND=ireals), ALLOCATABLE :: a_p_w2 (:,:)
  REAL (KIND=ireals), ALLOCATABLE :: a_w_T  (:,:)
  REAL (KIND=ireals), ALLOCATABLE :: a_w_w1 (:,:)
  REAL (KIND=ireals), ALLOCATABLE :: a_w_w2 (:,:)
  REAL (KIND=ireals), ALLOCATABLE :: a_T_w1 (:,:)
  REAL (KIND=ireals), ALLOCATABLE :: a_T_w2 (:,:)

  REAL (KIND=ireals) :: acrlatg
  REAL (KIND=ireals) :: one_m_wgtfac
  REAL (KIND=ireals) :: dts_inv

  INTEGER            :: i, j, k
  INTEGER            :: izstata

  ALLOCATE( a_w_p1      (1:ie,2:ke), STAT=izstata )
  ALLOCATE( a_w_p2      (1:ie,2:ke), STAT=izstata )
  ALLOCATE( a_w_T       (1:ie,2:ke), STAT=izstata )
  ALLOCATE( a_w_w1      (1:ie,2:ke), STAT=izstata )
  ALLOCATE( a_w_w2      (1:ie,1:ke), STAT=izstata )
  ALLOCATE( a_p_w1      (1:ie,1:ke), STAT=izstata )
  ALLOCATE( a_p_w2      (1:ie,1:ke), STAT=izstata )
  ALLOCATE( a_T_w1      (1:ie,1:ke), STAT=izstata )
  ALLOCATE( a_T_w2      (1:ie,1:ke), STAT=izstata )

  IF ( izdebug >= 10 ) THEN
    WRITE(*,*) "Subr. [lhs_of_tridiag_system_for_w] ..."
    WRITE(*,*) "    zl_3D_div_damping=", zl_3D_div_damping
  END IF

  dts_inv = 1.0_ireals / dts

  DO  j = jstart, jend

    ! coefficients for the implicit system:
    ! (a_y_x1  means that this coefficient comes from the 
    !  variable 'x' of the prognostic equation for 'y')

    DO k=2, ke
      DO  i = istart, iend

        a_w_p1 (i,k) = -rho_inv_k_avg(i,j,k) *sqrtg_r_w(i,j,k) *beta_s_2_impl
        a_w_p2 (i,k) =   g * beta_b_2_impl * p_inv_k_avg(i,j,k)
        a_w_T  (i,k) = - g * p0_pT0_k_avg(i,j,k) * beta_b_1_impl

      END DO
    END DO

    DO k=1, ke
      DO  i = istart, iend

        a_p_w1(i,k) = -gamma_p_tot(i,j,k) * sqrtg_r_s(i,j,k) * beta_s_4_impl
        a_p_w2(i,k) = -beta_b_3_impl * g * rho0(i,j,k)       * 0.5_ireals
        a_T_w1(i,k) = -R_by_cv *T_tot(i,j,k) *sqrtg_r_s(i,j,k) * beta_s_6_impl
        a_T_w2(i,k) = beta_b_4_impl * dT0dz(i,j,k)             * 0.5_ireals

      END DO
    END DO

    DO k=2, ke
      DO  i = istart, iend

        one_m_wgtfac = 1.0_ireals-wgtfac(i,j,k)

        lgs_A(i,j,k) =                                        &
              a_w_p1(i,k) * ( -a_p_w1(i,k-1)                  &
                             + a_p_w2(i,k-1) )                &
             + one_m_wgtfac *                                 &
                   (   a_w_p2(i,k) * ( a_p_w1(i,k-1)          &
                                     - a_p_w2(i,k-1)   )      &
                     + a_w_T(i,k)  * ( a_T_w1(i,k-1)          &
                                     - a_T_w2(i,k-1) ) )

        lgs_B(i,j,k) = dts_inv**2                                 &
             + a_w_p1(i,k) * ( ( a_p_w1(i,k) + a_p_w1(i,k-1) )    &
                             - ( a_p_w2(i,k) - a_p_w2(i,k-1) ) )  &
             + wgtfac(i,j,k) *                                    &
                   (   a_w_p2(i,k) * ( a_p_w1(i,k)                &
                                     - a_p_w2(i,k) )              &
                     + a_w_T(i,k)  * ( a_T_w1(i,k)                &
                                     - a_T_w2(i,k)   ) )          &
             - one_m_wgtfac *                                     &
                   (   a_w_p2(i,k) * ( a_p_w1(i,k-1)              &
                                     + a_p_w2(i,k-1) )            &
                     + a_w_T(i,k)  * ( a_T_w1(i,k-1)              &
                                     + a_T_w2(i,k-1) ) )

        lgs_C(i,j,k) =                                 &
             - a_w_p1(i,k) * ( a_p_w1(i,k)             &
                             + a_p_w2(i,k) )           &
             - wgtfac(i,j,k) *                         &
                    ( a_w_p2(i,k) * ( a_p_w1(i,k)      &
                                    + a_p_w2(i,k) )    &
                    + a_w_T(i,k)  * ( a_T_w1(i,k)      &
                                    + a_T_w2(i,k) ) )

      END DO
    END DO

    IF ( l_3D_div_damping ) THEN

      DO k=2, ke
        DO  i = istart, iend
          a_w_w1 (i,k) =  rho_inv_k_avg(i,j,k) * ( -sqrtg_r_w(i,j,k) )   
        END DO
      END DO

      DO k=1, ke
        DO  i = istart, iend
          a_w_w2 (i,k) = alpha_div_v_to_h * alpha_div_h(i,j,k) * rho(i,j,k)  &
                            * sqrtg_r_s(i,j,k) * beta_d_4_impl
        END DO
      END DO

      DO k=2, ke
        DO  i = istart, iend

          lgs_A(i,j,k) = lgs_A(i,j,k)                    &
               + dts_inv * a_w_w1(i,k) * a_w_w2(i,k-1)

          lgs_B(i,j,k) = lgs_B(i,j,k)                    &
               - dts_inv * a_w_w1(i,k) * ( a_w_w2(i,k) + a_w_w2(i,k-1) )

          lgs_C(i,j,k) = lgs_C(i,j,k)                    &
               + dts_inv * a_w_w1(i,k) * a_w_w2(i,k)

        END DO
      END DO

    END IF

  END DO

  ! -- Boundary conditions for w:

  DO  j = jstart, jend
    DO  i = istart, iend

      ! at top boundary: (w=0)
      lgs_A  (i,j,1) = 0.0_ireals
      lgs_B  (i,j,1) = 1.0_ireals
      lgs_C  (i,j,1) = 0.0_ireals

      ! at bottom boundary (d zeta/dt=0)
      lgs_A  (i,j,ke1) = 0.0_ireals
      lgs_B  (i,j,ke1) = 1.0_ireals
      lgs_C  (i,j,ke1) = 0.0_ireals

    END DO
  END DO


  IF ( .NOT. l_use_tridag ) THEN
    ! 1st step in LU decomposition:
    ! transform into a bidiagonal system:
    !   a'(k) w(k-1) + 1 * w(k) + 0 * w(k+1) = d'(k)
    ! (remark: lgs_B is used as an auxiliary variable)

    DO  j = jstart, jend
      DO  i = istart, iend
        lgs_B(i,j,ke) = 1.0_ireals / lgs_B(i,j,ke)
        lgs_A(i,j,ke) = lgs_A(i,j,ke) * lgs_B(i,j,ke)      
      END DO
    END DO

    DO k = ke-1, 2, -1
      DO  j = jstart, jend
        DO  i = istart, iend
          lgs_B(i,j,k) = 1.0_ireals/(lgs_B(i,j,k)-lgs_C(i,j,k)*lgs_A(i,j,k+1))
          lgs_A(i,j,k) = lgs_A(i,j,k) * lgs_B(i,j,k)      
        END DO
      END DO
    END DO

  ELSE
    ! nothing to do here
  END IF


  DEALLOCATE( a_w_p1  )
  DEALLOCATE( a_w_p2 )
  DEALLOCATE( a_w_T )
  DEALLOCATE( a_w_w1  )
  DEALLOCATE( a_w_w2 )
  DEALLOCATE( a_p_w1  )
  DEALLOCATE( a_p_w2  )
  DEALLOCATE( a_T_w1  )
  DEALLOCATE( a_T_w2  )

END SUBROUTINE lhs_of_tridiag_system_for_w

!=============================================================================
!=============================================================================

SUBROUTINE bottom_BC_for_w( u_sfc, v_sfc, dts, lgs_rhs )

  !--------------------------------------------------------------------------
  ! Description:
  !
  !   free slip bottom boundary condition for w:  d zeta / dt = 0 
  !
  ! Input:
  !   u_sfc(:,:), v_sfc(:,:)
  !   dts
  !   hhl
  !
  ! Output:
  !   lgs_rhs(:,:,ke1)
  !      (remark: the appropriate lgs_A, lgs_B, lgs_C(:,:,ke1) are set 
  !       in subr. 'lhs_of_tridiag_system_for_w')
  !
  ! Methods:
  !   according to itype_bbc_w
  !    =0/1: RK-like method following iadv_order (not yet implemented)
  !    =2/3: differencing following iadv_order without RK stepping
  !    =4/5: 4th-order centered differences
  !    =6/7: 2nd-order centered differences
  !    =0/2/4/6: include extrapolation of horizontal wind to sfc
  !    =1/3/5/7: no extrapolation of horizontal wind to sfc
  !  or new nomenclature in the form "1ed"
  !
  !--------------------------------------------------------------------------

  USE data_runcontrol, ONLY:  &
    iadv_order    ! order of the advection scheme in dyn. core

  USE src_advection_rk, ONLY:  &
    horiz_adv_driver

  IMPLICIT NONE

  REAL (KIND=ireals), INTENT(in)  :: u_sfc(ie,je), v_sfc(ie,je)
  REAL (kind=ireals), INTENT(in)  :: dts

  REAL (KIND=ireals), INTENT(out) :: lgs_rhs(ie,je)

  INTEGER            :: i, j
  REAL (KIND=ireals), ALLOCATABLE :: dlam_dt(:,:), dphi_dt(:,:)
  REAL (KIND=ireals), ALLOCATABLE :: tend(:,:)
  INTEGER            :: istat
  REAL (KIND=ireals) :: ssign
  REAL (KIND=ireals) :: r_earth_recip

  CHARACTER(LEN=100) :: yzerrmsg
  INTEGER            :: izerror 

  ALLOCATE( dlam_dt(ie,je), STAT=istat )
  ALLOCATE( dphi_dt(ie,je), STAT=istat )
  ALLOCATE( tend   (ie,je), STAT=istat )

  ssign = SIGN( 1.0_ireals, dts )

  r_earth_recip = 1.0_ireals / r_earth

  !------------------------------------------------------------------
  ! step 1:  at first, determine values u_sfc/v_sfc
  !   for u at the uw-position at the surface (at (i+1/2,j,ke+1/2) ) and
  !   for v at the vw-position at the surface (at (i,j+1/2,ke+1/2) ) :
  !------------------------------------------------------------------

  ! CALL calc_uv_surface( u, v )

  !--------------------------------------------------------------
  ! step 2: set lgs-coefficients for the free slip bottom BC
  !--------------------------------------------------------------

  ! set lhs of the equation system for w:
  ! (is done already in subr. 'lhs_of_tridiag_system_for_w')
  !lgs_A  (:,:,ke1) = 0.0_ireals
  !lgs_B  (:,:,ke1) = 1.0_ireals
  !lgs_C  (:,:,ke1) = 0.0_ireals

  ! set rhs of the equation system for w:
  IF ( (itype_bbc_w==2) .OR. (itype_bbc_w==3)  ) THEN

    ! calculate the rhs of w = ... via the advection operator for hhl

    DO  j = jstart, jend
      DO  i = istart, iend
        dlam_dt(i,j) = 0.5_ireals * ( u_sfc(i,j) + u_sfc(i-1,j) )    &
                           * acrlat(j,1)
        dphi_dt(i,j) = 0.5_ireals * ( v_sfc(i,j) + v_sfc(i,j-1) )    &
                           * r_earth_recip
      END DO
    END DO

    tend(:,:) = 0.0_ireals

    CALL horiz_adv_driver( hhl(:,:,ke1), dlam_dt, dphi_dt, tend, ssign,    &
      istart, iend, jstart, jend, 1, 1, 1, 1, iadv_order )

    lgs_rhs(:,:) = - tend(:,:)


  ELSE IF ( (itype_bbc_w==4) .OR. (itype_bbc_w==5) &
       .OR. (i_bbc_w_dhdx_order == 4) ) THEN

    ! hhl with centered differences 4th order
    DO  j = jstart, jend
      DO  i = istart, iend

        ! dz_dlam_4ord, ... are at the w-positions:
        lgs_rhs(i,j) = 0.5_ireals                                            &
            * (dz_dlam_4ord(i,j) *(u_sfc(i,j) + u_sfc(i-1,j)) *acrlat(j,1)   &
              +dz_dphi_4ord(i,j) *(v_sfc(i,j) + v_sfc(i,j-1)) *r_earth_recip )
      END DO
    END DO

  ELSE IF ( (itype_bbc_w==6) .OR. (itype_bbc_w==7) &
       .OR. (i_bbc_w_dhdx_order == 2) ) THEN

    ! hhl with centered differences 2nd order
    DO  j = jstart, jend
      DO  i = istart, iend
        !lgs_rhs(i,j) = 0.5_ireals                                         &
        !     * ( ( dz_dlam(i,  j,  ke1) * u_sfc(i  ,j  )                  &
        !         + dz_dlam(i-1,j,  ke1) * u_sfc(i-1,j  ) ) * acrlat(j,1)  &
        !       + ( dz_dphi(i,  j,  ke1) * v_sfc(i,  j  )                  &
        !         + dz_dphi(i,  j-1,ke1) * v_sfc(i,  j-1) ) * r_earth_recip )

        ! this is more consistent with the definition of dzeta/dt above:
        lgs_rhs(i,j) = acrlat(j,1) * 0.5_ireals                             &
             * ( ( dz_dlam(i,  j,  ke1) * u_sfc(i  ,j  )                    &
                 + dz_dlam(i-1,j,  ke1) * u_sfc(i-1,j  ) )                  &
               + ( dz_dphi(i,  j,  ke1) * v_sfc(i,  j  ) * crlat(j  ,2)     &
                 + dz_dphi(i,  j-1,ke1) * v_sfc(i,  j-1) * crlat(j-1,2) ) )
      END DO
    END DO

  ELSE IF ( (itype_bbc_w==8) .OR. (itype_bbc_w==9) &
       .OR. (i_bbc_w_dhdx_order == 3) ) THEN

    ! calculate the rhs of w = .. via the 3rd order advection operator for hhl

    DO  j = jstart, jend
      DO  i = istart, iend
        dlam_dt(i,j) = 0.5_ireals * ( u_sfc(i,j) + u_sfc(i-1,j) )    &
                           * acrlat(j,1)
        dphi_dt(i,j) = 0.5_ireals * ( v_sfc(i,j) + v_sfc(i,j-1) )    &
                           * r_earth_recip
      END DO
    END DO

    tend(:,:) = 0.0_ireals

    CALL horiz_adv_driver( hhl(:,:,ke1), dlam_dt, dphi_dt, tend, ssign,    &
      istart, iend, jstart, jend, 1, 1, 1, 1, 3 )

    lgs_rhs(:,:) = - tend(:,:)

  ELSE
    yzerrmsg="(step 2) this case for itype_bbc_w is not available!"
    izerror=100
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'bottom_BC_for_w')
  END IF


  DEALLOCATE( dlam_dt )
  DEALLOCATE( dphi_dt )
  DEALLOCATE( tend  )

END SUBROUTINE bottom_BC_for_w

!=============================================================================
!=============================================================================

SUBROUTINE tridag( a, b, c, rhs, y, n)

  !--------------------------------------------
  !
  ! description:
  ! solve the tridiagonal linear system of equations
  !    a(k) y(k-1) + b(k) y(k) + c(k) y(k+1) = rhs(k), k=1, ..., n
  ! (--> a(1)=0 and c(n)=0 assumed (i.e. they are not used)
  !
  ! input 
  !   a, b, c : 3 diagonals of the matrix (are not modified)
  !   rhs : (is not modified)
  !   n : number of equations
  !
  ! method:
  ! Fortran-version of Numerical recipes in C, 2nd edition
  ! (optimized for vector computers, should be fast on scalar machines, too)
  !
  !-------------------------------------------------------

  USE data_modelconfig, ONLY:  &
    ie, je,                    &
    istart, iend,              &
    jstart, jend

  IMPLICIT NONE

  INTEGER (KIND=iintegers), INTENT(IN) :: n

  REAL (KIND=ireals), INTENT(IN) :: a  (1:ie, 1:je, 1:n)
  REAL (KIND=ireals), INTENT(IN) :: b  (1:ie, 1:je, 1:n)
  REAL (KIND=ireals), INTENT(IN) :: c  (1:ie, 1:je, 1:n)
  REAL (KIND=ireals), INTENT(IN) :: rhs(1:ie, 1:je, 1:n)

  REAL (KIND=ireals), INTENT(OUT) :: y(1:ie, 1:je, 1:n)

  REAL (KIND=ireals), ALLOCATABLE :: work(:,:,:)
  REAL (KIND=ireals), ALLOCATABLE :: bet(:,:)

  INTEGER (KIND=iintegers) :: i,j,k
  INTEGER (KIND=iintegers) :: istat

  ALLOCATE( work(1:ie, 1:je, 1:n ), STAT=istat )
  ALLOCATE( bet (1:ie, 1:je),       STAT=istat )

  !WRITE(*,*) "Subr. [tridag] ..."

  ! do i=1, ie
  !   do j=1, je
  !     if ( b(i,j,1) == 0.0_ireals ) then
  !       write(*,*) "Error 1 in tridag"
  !     end if
  !   end do
  ! end do

  DO j=jstart, jend
    DO i=istart, iend
      bet(i,j) = b(i,j,1)
      y(i,j,1) = rhs(i,j,1) / bet(i,j)
    END DO
  END DO

  DO k=2, n
    DO j=jstart, jend
      DO i=istart, iend
        work(i,j,k) = c(i,j,k-1) / bet(i,j)
        bet(i,j) = b(i,j,k) - a(i,j,k) * work(i,j,k)
        y(i,j,k) = ( rhs(i,j,k) - a(i,j,k) * y(i,j,k-1) )/bet(i,j)
      END DO
    END DO
  END DO

  ! backsubstitution:
  DO k=n-1, 1, -1
    DO j=jstart, jend
      DO i=istart, iend
        y(i,j,k) = y(i,j,k) - work(i,j,k+1) * y(i,j,k+1)
      END DO
    END DO
  END DO

  DEALLOCATE( work )
  DEALLOCATE( bet)

END SUBROUTINE tridag

!=============================================================================
!=============================================================================


SUBROUTINE init_div_damping_coeff( xkd, dts, alpha_div_h, alpha_div_v_to_h )

  !---------------------------------------------------------------------------
  !
  ! Description:
  ! calculate the horizontal divergence damping coefficient according to 
  ! the dim.less value xkd and the slope of the orography
  !
  !---------------------------------------------------------------------------

  USE data_modelconfig, ONLY :   &
    dlon, dlat
  USE data_constants,   ONLY :   &
    pi

  IMPLICIT NONE

  REAL (Kind=ireals), INTENT(in)  :: xkd, dts
  REAL (Kind=ireals), INTENT(out) :: alpha_div_h(ie,je,ke)
  REAL (Kind=ireals), INTENT(out) :: alpha_div_v_to_h

  REAL (Kind=ireals)  :: alpha_div_h_limit
  REAL (Kind=ireals)  :: alpha_div_h_min, alpha_div_h_max
  REAL (Kind=ireals)  :: alpha_div_h_avg
  REAL (Kind=ireals)  :: c_sound = 330.0_ireals
  REAL (kind=ireals)  :: work(ie,je) 

  REAL (Kind=ireals)  :: dlon_arg, dlat_arg
  REAL (Kind=ireals)  :: delta_x_inv, delta_x_inv_0, delta_y_inv
  REAL (Kind=ireals)  :: ax, ay, h1

  REAL (Kind=ireals)  :: dz_dlam, dz_dphi
  REAL (Kind=ireals)  :: wght0, wght1, wght2, wght_fac

  INTEGER :: i,j,k

  INTEGER (KIND=iintegers) ::  &
    kzdims(24)          ! Vertical dimensions for exchg_boundaries
  INTEGER :: izerror

  CHARACTER (LEN=100) :: yzerrmsg


  wght_fac = 0.25_ireals
  wght2 = wght_fac**2
  wght1 = wght_fac * ( 1.0_ireals - 2.0_ireals * wght_fac )
  wght0 = ( 1.0_ireals - 2.0_ireals * wght_fac )**2

  alpha_div_h_limit =  xkd * c_sound**2 * dts

  dlon_arg = dlon * pi / 180.0_ireals
  dlat_arg = dlat * pi / 180.0_ireals

  delta_x_inv_0 = 1.0_ireals / ( r_earth * dlon_arg )
  delta_y_inv   = 1.0_ireals / ( r_earth * dlat_arg )

  alpha_div_h_max = -1.0e30
  alpha_div_h_min =  1.0e30

  alpha_div_h_avg = 0.0_ireals

  IF ( .NOT. limit_alpha_div_by_slope_stability ) THEN

    alpha_div_h(:,:,:) = alpha_div_h_limit

  ELSE   ! limit_alpha_div_by_slope_stability = .TRUE.

    alpha_div_h(:,:, 1:vcoord%kflat-1) = alpha_div_h_limit

    DO k=vcoord%kflat, ke

      DO j=1, je-1
        delta_x_inv = delta_x_inv_0 / crlat(j,1)
        DO i=1, ie-1

          ! fields dz_dlam, dz_dphi are not yet available here:
          dz_dlam = ( hhl(i+1,j,k) - hhl(i,j,k) ) * eddlon
          dz_dphi = ( hhl(i,j+1,k) - hhl(i,j,k) ) * eddlat
          ax = (delta_x_inv *(2.0_ireals + dlon_arg *ABS(dz_dlam) ) )**2
          ay = (delta_y_inv *(2.0_ireals + dlat_arg *ABS(dz_dphi) ) )**2

          !ax = (delta_x_inv *(2.0_ireals + dlon_arg *ABS(dz_dlam(i,j,k))))**2
          !ay = (delta_y_inv *(2.0_ireals + dlat_arg *ABS(dz_dphi(i,j,k))))**2

          !work(i,j) = 5.0_ireals * 2.0_ireals / ( dts * ( ax + ay ) )
          work(i,j) = divdamp_slope * 2.0_ireals / ( dts * ( ax + ay ) )

        END DO
      END DO

      DO j=2, je-2
        DO i=2, ie-2

          ! horizontal averaging (over 9 points):
          h1 =   wght2 * ( work(i+1,j+1) + work(i+1,j-1)     &
                         + work(i-1,j+1) + work(i-1,j-1) )   &
               + wght1 * ( work(i+1,j)   + work(i-1,j)       &
                         + work(i  ,j+1) + work(i  ,j-1) )   &
               + wght0 *   work(i,j)

          alpha_div_h(i,j,k) = MIN( h1, alpha_div_h_limit)

        END DO
      END DO

      alpha_div_h(1   ,:,k) = alpha_div_h(2   ,:,k)
      alpha_div_h(ie-1,:,k) = alpha_div_h(ie-2,:,k)
      alpha_div_h(ie  ,:,k) = alpha_div_h(ie-2,:,k)
      alpha_div_h(:   ,1,k) = alpha_div_h(:,   2,k)
      alpha_div_h(:,je-1,k) = alpha_div_h(:,je-2,k)
      alpha_div_h(:,je  ,k) = alpha_div_h(:,je-2,k)

    END DO

  END IF

  ! ----- Exchange of alpha_div ---------------

  kzdims(1:24) =                          &
       (/  ke,   0,   0,   0,   0,      &
            0,   0,   0,   0,   0,      &
            0,   0,   0,   0,   0,      &
            0,   0,   0,   0,   0,      &
            0,   0,   0,   0  /)

  ! non-optimized version (icase=0 <--> ldatatypes=.FALSE.)
  CALL exchg_boundaries                                                      &
      (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,  &
       kzdims, jstartpar, jendpar, nboundlines, nboundlines, my_cart_neigh,  &
       lperi_x, lperi_y, l2dim,                                              &
       25000+nexch_tag, .FALSE., ncomm_type, izerror, yzerrmsg,              &
       alpha_div_h(:,:,:) )

  IF ( izerror /= 0_iintegers ) THEN
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'fast_waves_strong_conserv')
  END IF


  ! choose ratio between horizontal and vertical divergence damping coeff.:
  IF ( l_3D_div_damping ) THEN 
    alpha_div_v_to_h = 1.0_ireals
  ELSE
    alpha_div_v_to_h = 0.0_ireals
  END IF

END SUBROUTINE init_div_damping_coeff

!=============================================================================
!=============================================================================


SUBROUTINE calc_uv_surface( u, v, extrapol_uv_order, u_sfc, v_sfc )

  !--------------------------------------------------------------------------
  ! Description:
  !
  !   calculate values for u and v at the bottom surface 
  !   (i.e. in (i+1/2, j, ke+1/2), (i, j+1/2, k+1/2), respectively)
  !
  ! Input:
  !   u(:,:,:), v(:,:,:)
  !   extrapol_uv_order
  !
  ! Output:
  !   u_sfc(:,:), v_sfc(:,:)
  !
  ! Method:
  !   use 0th, 1st or 2nd order vertical extrapolation dependent on 
  !   the variable 'extrapol_uv_order'
  !--------------------------------------------------------------------------

  USE data_runcontrol, ONLY:  &
    ltur            ! if .TRUE., turbulence parameterization is used
  !lfreeslip_sfc   ! if .TRUE., surface friction is turned off despite using
  !                ! a turbulence scheme (relevant for idealized simulations)

  USE src_artifdata, ONLY:  &
    lnosurffluxes_m  ! ersetzt lfreeslip_sfc

  USE grid_metrics_utilities, ONLY:  &
    wgtfac_u, wgtfac_v, wgtfacq_u, wgtfacq_v

  REAL (KIND=ireals), INTENT(in)  :: u(ie,je,ke)
  REAL (KIND=ireals), INTENT(in)  :: v(ie,je,ke)
  INTEGER,            INTENT(in)  :: extrapol_uv_order

  REAL (KIND=ireals), INTENT(out) :: u_sfc(ie,je)
  REAL (KIND=ireals), INTENT(out) :: v_sfc(ie,je)

  INTEGER :: i, j

  CHARACTER(LEN=80) :: yzerrmsg
  INTEGER :: izerror

  IF ( extrapol_uv_order == 0 ) THEN
    ! free slip condition (no extrapolation):
    u_sfc(:,:) = u(:,:,ke)
    v_sfc(:,:) = v(:,:,ke)

  ELSE IF ( extrapol_uv_order == 1 ) THEN
    ! Linear extrapolation for bottom level if u/v(ke) is used
    ! for bottom boundary condition of w (by G. Zaengl)
    DO j = jstart-1, jend+1
      !CDIR ON_ADB(u)
      !CDIR ON_ADB(v)
      DO i = istart-1, iend+1
        u_sfc(i,j) = u(i,j,ke) + wgtfac_u(i,j,ke1) &
                                    * (u(i,j,ke) - u(i,j,ke-1))
        v_sfc(i,j) = v(i,j,ke) + wgtfac_v(i,j,ke1) &
                                    * (v(i,j,ke) - v(i,j,ke-1))
      ENDDO
    ENDDO

  ELSE IF ( extrapol_uv_order == 2 ) THEN
    ! Quadratic extrapolation (by G. Zaengl):
    DO j = jstart-1, jend+1
      !CDIR ON_ADB(u)
      !CDIR ON_ADB(v)
      DO i = istart-1, iend+1
        u_sfc(i,j) = u(i,j,ke  ) * wgtfacq_u(i,j,1)       &
                   + u(i,j,ke-1) * wgtfacq_u(i,j,2)       &
                   + u(i,j,ke-2) * wgtfacq_u(i,j,3)
        v_sfc(i,j) = v(i,j,ke  ) * wgtfacq_v(i,j,1)       &
                   + v(i,j,ke-1) * wgtfacq_v(i,j,2)       &
                   + v(i,j,ke-2) * wgtfacq_v(i,j,3)

        ! Limit extrapolated value to the absolute value at the lowermost 
        ! main level; this is needed to avoid spurious overshoots 
        ! for underresolved katabatic flows
        IF ( ltur .AND. .NOT. lnosurffluxes_m ) THEN
          IF ( u(i,j,ke) > 0.0_ireals ) THEN
            u_sfc(i,j) = MIN( u_sfc(i,j), u(i,j,ke) )
          ELSE
            u_sfc(i,j) = MAX( u_sfc(i,j), u(i,j,ke) )
          ENDIF
          IF ( v(i,j,ke) > 0.0_ireals ) THEN
            v_sfc(i,j) = MIN( v_sfc(i,j), v(i,j,ke) )
          ELSE
            v_sfc(i,j) = MAX( v_sfc(i,j), v(i,j,ke) )
          ENDIF
        ENDIF
        ! Ensure that extrapolation does not overshoot to opposite sign
        ! In that case, surface wind is set to zero
        IF ( u_sfc(i,j)*u(i,j,ke) <= 0.0_ireals) &
          u_sfc(i,j) = 0.0_ireals
        IF ( v_sfc(i,j)*v(i,j,ke) <= 0.0_ireals) &
          v_sfc(i,j) = 0.0_ireals
      ENDDO
    ENDDO

  ELSE
    yzerrmsg="this case for extrapol_uv_order is not available!"
    izerror=100
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'calc_uv_surface')
  END IF

END SUBROUTINE calc_uv_surface

!=============================================================================
!=============================================================================

SUBROUTINE calc_Z_horiz ( u, v, u_sfc, v_sfc, Z_horiz )

  !------------------------------------------------------
  ! Description:
  !   calculate (at w-position =(i, j, k-1/2) ):
  !   ( the 'horizontal' contribution of
  !       Z = r * cos phi * sqrt(g) * dzeta/dt )
  !   Z_horiz =       Z_x    +      Z_y
  !           =  dz/dlam * u + dz/dphi * cos phi * v 
  !   or more precisely:
  !           =   avg_lam ( dz/dlam * avg_z[u] )
  !             + avg_phi ( dz/dphi * cos(phi) * avg_z[v] )
  !
  ! Input:
  !   u, v
  !   u_sfc, v_sfc
  !
  ! Output:
  !   Z_horiz
  !------------------------------------------------------

  USE grid_metrics_utilities, ONLY:  &
    wgtfac_u, wgtfac_v

  IMPLICIT NONE

  REAL (KIND=ireals), INTENT (IN) ::    &
    u (ie,je,ke),   &
    v (ie,je,ke)

  REAL (KIND=ireals), INTENT (IN) ::    &
    u_sfc (ie,je),   &
    v_sfc (ie,je)

  REAL (KIND=ireals), INTENT (OUT) ::    &
    Z_horiz(ie,je,ke1)

  REAL (KIND=ireals), ALLOCATABLE :: Z_x(:,:), Z_y(:,:)

  INTEGER :: i, j, k
  INTEGER :: istat

  REAL (KIND=ireals) :: r_earth_recip

  CHARACTER(LEN=100) yzerrmsg

  ALLOCATE( Z_x (ie, je), STAT=istat)
  ALLOCATE( Z_y (ie, je), STAT=istat)

  r_earth_recip = 1.0 / r_earth

  Z_horiz(:,:,1:MAX(1,vcoord%kflat)) = 0.0_ireals

  DO k=MAX(2,vcoord%kflat+1), ke

    DO j=jstart, jend
      DO i=istart-1, iend
        Z_x(i,j) =   dz_dlam(i,j,k) *                          &
             (                wgtfac_u(i,j,k)   * u(i,j,k  )   &
             + ( 1.0_ireals - wgtfac_u(i,j,k) ) * u(i,j,k-1) )
      END DO
    END DO

    DO j=jstart-1, jend
      DO i=istart, iend
        Z_y(i,j) =   dz_dphi(i,j,k) * crlat(j,2) *             &
             (                wgtfac_v(i,j,k)   * v(i,j,k  )   &
             + ( 1.0_ireals - wgtfac_v(i,j,k) ) * v(i,j,k-1) )
      END DO
    END DO

    DO j=jstart, jend
      DO i=istart, iend
        Z_horiz(i,j,k) =                                   & 
            0.5_ireals * ( ( Z_x(i,j) + Z_x(i-1,j  ) )     &
                         + ( Z_y(i,j) + Z_y(i  ,j-1) ) )
      END DO
    END DO

  END DO

  ! --- Boundary treatment: ------

  ! top boundary (k=1) is flat in any case (see above)

  ! at bottom boundary (k=ke+1): use extrapolated values for u, v

  IF ( itype_bound_for_Z_horiz == 2 ) THEN

    k=ke1

    DO j=jstart, jend
      DO i=istart-1, iend
        Z_x(i,j) = dz_dlam(i,j,ke1) * u_sfc(i,j)
      END DO
    END DO

    DO j=jstart-1, jend
      DO i=istart, iend
        Z_y(i,j) = dz_dphi(i,j,ke1) * crlat(j,2) * v_sfc(i,j)
      END DO
    END DO

    DO j=jstart, jend
      DO i=istart, iend
        Z_horiz(i,j,ke1) = 0.5_ireals            &
             * ( ( Z_x(i,j) + Z_x(i-1,j  ) )     &
               + ( Z_y(i,j) + Z_y(i  ,j-1) ) )
      END DO
    END DO

  ELSE IF ( itype_bound_for_Z_horiz == 4 ) THEN

    k=ke1

    DO j=jstart, jend
      DO i=istart, iend
        Z_horiz(i,j,ke1) = 0.5_ireals                                        &
            * ( dz_dlam_4ord(i,j) *(u_sfc(i,j) + u_sfc(i-1,j))               &
              + dz_dphi_4ord(i,j) *(v_sfc(i,j) + v_sfc(i,j-1)) * crlat(j,1) )

        ! alternative:
        !Z_horiz(i,j,ke1) = 0.5_ireals                                      &
        !     * ( dz_dlam_4ord(i,j) * ( u_sfc(i,j) + u_sfc(i-1,j) )         &
        !       + dz_dphi_4ord(i,j) * ( v_sfc(i,j)   * crlat(j,  2)         &
        !                             + v_sfc(i,j-1) * crlat(j-1,2) ) )
      END DO
    END DO

  ELSE
    yzerrmsg = "false value in itype_bound_for_Z_horiz"
    CALL model_abort (my_cart_id, 1234, yzerrmsg, "calc_Z_horiz")
  END IF


  ! at top boundary (k=1): free-slip condition for dzeta/dt
  ! here: no BC is set
  !Z_horiz(:,:,1)   = 1.0e20_ireals

  ! at bottom boundary (k=ke+1): explicit (!) free-slip condition for dzeta/dt
  ! here: no BC is set
  !Z_horiz(:,:,ke1) = 1.0e20_ireals


  IF ( l2dim ) THEN
    DO j=1, nboundlines
      Z_horiz(:,jstart-j,:) = Z_horiz(:,jstart,:)
      Z_horiz(:,jend  +j,:) = Z_horiz(:,jstart,:)
    END DO
  END IF


  DEALLOCATE( Z_x )
  DEALLOCATE( Z_y )

END SUBROUTINE calc_Z_horiz

!=============================================================================
!=============================================================================

SUBROUTINE calrho_fw ( T_tot, gamma_p_tot, qv, q_cond, rho )

  !---------------------------------------------------------------------------
  ! Description:
  !   This routine computes the air density for present time level.
  !
  ! Method:
  !   Application of the gas-law, rho = p/(r_d*tv)  
  !
  !---------------------------------------------------------------------------
  !

  USE data_constants  , ONLY :   &
    R_d,          & ! gas constant for dry air
    rvd_m_o,      & ! r_v/r_d - 1; hier fuer Subr. calrho
    gamma


  REAL (KIND=ireals), INTENT (IN)          ::    &
    T_tot       (ie,je,ke),   & ! temperature
    gamma_p_tot (ie,je,ke),   & ! = c_p / c_v * total pressure
    qv          (ie,je,ke),   & ! specific humidity
    q_cond      (ie,je,ke)      ! specific cloud water content

  REAL (KIND=ireals), INTENT (OUT)       ::    &
    rho  (ie,je,ke)      ! air density 

  REAL (KIND=ireals)  :: gamma_Rd    ! = c_p / c_v * R_d

  !---------------------------------------------------------------------------
  ! Begin subroutine calrho

  gamma_Rd = gamma * r_d

  rho(:,:,:) = gamma_p_tot(:,:,:) / ( gamma_Rd * T_tot(:,:,:) *              &
                     ( 1.0_ireals + rvd_m_o * qv(:,:,:) - q_cond(:,:,:) ) )

END SUBROUTINE calrho_fw



REAL (KIND=ireals) FUNCTION rho_fw_gp ( T_tot, p_tot, qv, q_cond )

  !---------------------------------------------------------------------------
  ! Description:
  !   This routine computes the air density for present time level.
  !   (similar to subr. 'calrho_fw_gp', but for a single gridpoint).
  !
  ! Method:
  !   Application of the gas-law, rho = p/(r_d*tv)  
  !
  !---------------------------------------------------------------------------
  !

  USE data_constants  , ONLY :   &
    R_d,          & ! gas constant for dry air
    rvd_m_o         ! r_v/r_d - 1; hier fuer Subr. calrho

  REAL (KIND=ireals), INTENT (IN)          ::    &
    T_tot       ,   & ! temperature
    p_tot       ,   & ! total pressure
    qv          ,   & ! specific humidity
    q_cond            ! specific cloud water content

  rho_fw_gp = p_tot / ( R_d * T_tot * ( 1.0_ireals + rvd_m_o * qv - q_cond ) )

END FUNCTION rho_fw_gp

!=============================================================================
!=============================================================================

SUBROUTINE radiative_lbc( u, v, w, pp, Tp, dts,                   &
  bd_west, bd_east, bd_north, bd_south )

  !--------------------------------------------------------------------------
  ! Description:
  !   Radiative lateral boundary condition for all dynamical variables
  !
  ! Input variables:
  !   dynamical fields 
  !     u(:,:,:), v(:,:,:), w(:,:,:), pp(:,:,:), Tp(:,:,:)
  !
  ! Output variables:
  !   2-dim. boundary fields for the 5 variables 
  !      bd_west(:,:,5),  bd_east(:,:,5),  bd_north(:,:,5),  bd_south(:,:,5) 
  !
  ! Method:
  !   solution of a one-way wave equation
  !   df/dt + c * df/dx = 0     at the boundary
  !       c ~ sound velocity
  !--------------------------------------------------------------------------

  REAL (kind=ireals), INTENT(in)  ::  u(ie,je,ke)
  REAL (kind=ireals), INTENT(in)  ::  v(ie,je,ke)
  REAL (kind=ireals), INTENT(in)  ::  w(ie,je,ke1)
  REAL (kind=ireals), INTENT(in)  :: pp(ie,je,ke)
  REAL (kind=ireals), INTENT(in)  :: Tp(ie,je,ke)

  REAL (kind=ireals), INTENT(in)  :: dts

  REAL (kind=ireals), INTENT(out) :: bd_west ( je, ke1, 5)
  REAL (kind=ireals), INTENT(out) :: bd_east ( je, ke1, 5)
  REAL (kind=ireals), INTENT(out) :: bd_north( ie, ke1, 5)
  REAL (kind=ireals), INTENT(out) :: bd_south( ie, ke1, 5)

  REAL (kind=ireals), ALLOCATABLE :: cfl(:)  ! Courant-numbers
  REAL (kind=ireals) :: cs         ! sound velocity

  INTEGER :: i, j, k

  ! these values must be the same as in subr. 'set_radiative_lbc':
  INTEGER :: var_u  = 1
  INTEGER :: var_v  = 2
  INTEGER :: var_w  = 3
  INTEGER :: var_Tp = 4
  INTEGER :: var_pp = 5

  ! Indizes for the position of u, v, s (=scalars pp, Tp, or w).
  INTEGER i_west_u
  INTEGER i_west_v
  INTEGER i_west_s
  INTEGER i_east_u
  INTEGER i_east_v
  INTEGER i_east_s
  INTEGER j_north_u
  INTEGER j_north_v
  INTEGER j_north_s
  INTEGER j_south_u
  INTEGER j_south_v
  INTEGER j_south_s

  ALLOCATE( cfl(je) )

  WRITE(*,*) "subr. [radiative_lbc] ..."

  !cs = 300.0  ! rough estimate for the sound velocity
  cs = 50.0   ! rough estimate for the radiation velocity; 
  ! compromise between sound and gravity wave velocity

  !bd_west (:,:,:)=9999.0
  !bd_east (:,:,:)=9999.0
  !bd_north(:,:,:)=9999.0
  !bd_south(:,:,:)=9999.0

  ! Indizes for the boundary positions of u, v, s (=scalars pp, Tp, or w).
  ! They are the same as in subr. 'set_radiative_lbc':

  i_west_u = istartu-1
  i_west_v = istartv
  i_west_s = istart

  i_east_u = iendu+1
  i_east_v = iendv
  i_east_s = iend

  j_north_u = jendu
  j_north_v = jendv+1
  j_north_s = jend

  j_south_u = jstartu
  j_south_v = jstartv-1
  j_south_s = jstart


  ! western boundary
  IF ( my_cart_neigh(1) == -1 ) THEN

    IF ( l2dim ) THEN  
      cfl(:) = 1.0_ireals
    ELSE
      DO j=1, je
        cfl(j) = cs * dts * acrlat(j,1) * eddlon
      END DO
    END IF

    DO k=1, ke
      DO j = j_south_u, j_north_u
        ! Anm.: eigentlich nur:  j = j_south_u+1, j_north_u-1,
        ! aber so hat man auch den Fall  l2dim=.TRUE.  dabei
        ! (analog folgende Schleifen)

        bd_west(j,k, var_u) = u( i_west_u,   j,k)                       &
          &       + cfl(j) * (u( i_west_u+1, j,k) - u( i_west_u, j,k) )
      END DO
    END DO

    DO k=1, ke+1
      DO j = j_south_s, j_north_s
        bd_west(j,k, var_w) = w( i_west_s,   j,k)                       &
          &       + cfl(j) * (w( i_west_s+1, j,k) - w( i_west_s, j,k) )
      END DO
    END DO

    IF ( l2dim ) THEN
      DO k=1, ke
        DO j = j_south_s, j_north_s
          bd_west(j,k, var_pp) = pp( i_west_s,   j,k)                     &
            &        + cfl(j) * (pp( i_west_s+1, j,k) - pp( i_west_s, j,k) )
        END DO
      END DO
    ELSE
      DO k=1, ke
        DO j = j_south_s, j_north_s
          ! radiation velocity = 0 for pp!
          bd_west(j,k, var_pp) = pp( i_west_s, j,k)
        END DO
      END DO
    END IF

    DO k=1, ke
      DO j = j_south_s, j_north_s
        bd_west(j,k, var_Tp) = Tp( i_west_s,   j,k)                     &
          &        + cfl(j) * (Tp( i_west_s+1, j,k) - Tp( i_west_s, j,k) )
      END DO
    END DO

    IF ( .NOT. l2dim ) THEN
      DO j = 1, je
        cfl(j) = cs * dts * acrlat(j,2) * eddlon
      END DO

      DO k=1, ke
        DO j = j_south_v+1, j_north_v-1
          bd_west(j,k, var_v) = v( i_west_v,   j,k)                     &
            &       + cfl(j) * (v( i_west_v+1, j,k) - v( i_west_v, j,k) )
        END DO
      END DO
    ELSE
      bd_west(:,:, var_v) = 9999.0_ireals   ! senseless value
    END IF

  END IF


  ! eastern boundary
  IF ( my_cart_neigh(3) == -1 ) THEN

    IF ( l2dim ) THEN  
      cfl(:) = 1.0_ireals
    ELSE
      DO j=1, je
        cfl(j) = cs * dts * acrlat(j,1) * eddlon
      END DO
    END IF

    DO k=1, ke
      DO j = j_south_u, j_north_u
        ! Anm.: eigentlich nur:  j = j_south_u+1, j_north_u-1,
        ! aber so hat man auch den Fall  l2dim=.TRUE.  dabei
        ! (analog folgende Schleifen)

        bd_east(j,k, var_u) = u( i_east_u, j,k)                          &
          &       - cfl(j) * (u( i_east_u, j,k) - u( i_east_u-1, j,k) )
      END DO
    END DO

    DO k=1, ke1
      DO j = j_south_s, j_north_s
        bd_east(j,k, var_w) = w( i_east_s, j,k)                          &
          &       - cfl(j) * (w( i_east_s, j,k) - w( i_east_s-1, j,k) )
      END DO
    END DO

    IF ( l2dim ) THEN
      DO k=1, ke
        DO j = j_south_s, j_north_s
          bd_east(j,k, var_pp) = pp( i_east_s, j,k)                        &
            &        - cfl(j) * (pp( i_east_s, j,k) - pp( i_east_s-1, j,k) )
        END DO
      END DO
    ELSE
      DO k=1, ke
        DO j = j_south_s, j_north_s
          ! radiation velocity = 0 for pp!
          bd_east(j,k, var_pp) = pp( i_east_s, j,k)
        END DO
      END DO
    END IF

    DO k=1, ke
      DO j = j_south_s, j_north_s
        bd_east(j,k, var_Tp) = Tp( i_east_s, j,k)                        &
          &        - cfl(j) * (Tp( i_east_s, j,k) - Tp( i_east_s-1, j,k) )
      END DO
    END DO

    IF ( .NOT. l2dim ) THEN
      DO j=1, je
        cfl(j) = cs * dts * acrlat(j,2) * eddlon
      END DO

      DO k=1, ke
        DO j = j_south_v+1, j_north_v-1
          bd_east(j,k, var_v) = v( i_east_v, j,k)                        &
            &       - cfl(j) * (v( i_east_v, j,k) - v( i_east_v-1, j,k) )
        END DO
      END DO
    ELSE 
      bd_east(:,:, var_v) = 9999.0_ireals   ! senseless value
    END IF

  END IF


  ! northern boundary
  IF ( (my_cart_neigh(2) == -1) .AND. (.NOT. l2dim) ) THEN

    cfl(1) = cs * dts * eddlat / r_earth

    DO k=1, ke
      DO i = i_west_u+1, i_east_u-1
        bd_north(i,k, var_u) = u(i, j_north_u, k)                         &
          &        - cfl(1) * (u(i, j_north_u, k) - u(i, j_north_u-1, k) )
      END DO
    END DO

    DO k=1, ke
      DO i = i_west_v+1, i_east_v-1
        bd_north(i,k, var_v) = v(i, j_north_v, k)                         &
          &        - cfl(1) * (v(i, j_north_v, k) - v(i, j_north_v-1, k) )
      END DO
    END DO

    DO k=1, ke+1
      DO i = i_west_s+1, i_east_s-1
        bd_north(i,k, var_w) = w(i, j_north_s, k)                         &
          &        - cfl(1) * (w(i, j_north_s, k) - w(i, j_north_s-1, k) )
      END DO
    END DO

    IF ( l2dim ) THEN
      DO k=1, ke
        DO i = i_west_s+1, i_east_s-1
          bd_north(i,k, var_pp) = pp(i, j_north_s, k)                       &
            &         - cfl(1) * (pp(i, j_north_s, k) - pp(i, j_north_s-1, k) )
        END DO
      END DO
    ELSE
      DO k=1, ke
        DO i = i_west_s+1, i_east_s-1
          ! radiation velocity = 0 for pp!
          bd_north(i,k, var_pp) = pp(i, j_north_s, k)
        END DO
      END DO
    END IF

    DO k=1, ke
      DO i = i_west_s+1, i_east_s-1
        bd_north(i,k, var_Tp) = Tp(i, j_north_s, k)                        &
          &         - cfl(1) * (Tp(i, j_north_s, k) - Tp(i, j_north_s-1, k) )
      END DO
    END DO

  END IF


  ! southern boundary
  IF ( ( my_cart_neigh(4) == -1 ) .AND. (.NOT. l2dim) ) THEN

    cfl(1) = cs * dts * eddlat / r_earth

    DO k=1, ke
      DO i = i_west_u+1, i_east_u-1
        bd_south(i,k, var_u) = u(i, j_south_u,   k)                       &
          &        + cfl(1) * (u(i, j_south_u+1, k) - u(i, j_south_u, k) )
      END DO
    END DO

    DO k=1, ke
      DO i = i_west_v+1, i_east_v-1
        bd_south(i,k, var_v) = v(i, j_south_v, k)                         &
          &        + cfl(1) * (v(i, j_south_v+1, k) - v(i, j_south_v, k) )
      END DO
    END DO

    DO k=1, ke+1
      DO i = i_west_s+1, i_east_s-1
        bd_south(i,k, var_w) = w(i, j_south_s,   k)                       &
          &        + cfl(1) * (w(i, j_south_s+1, k) - w(i, j_south_s, k) )
      END DO
    END DO

    IF ( l2dim ) THEN
      DO k=1, ke
        DO i = i_west_s+1, i_east_s-1
          bd_south(i,k, var_pp) = pp(i, j_south_s,   k)                      &
            &         + cfl(1) * (pp(i, j_south_s+1, k) - pp(i, j_south_s, k) )
        END DO
      END DO
    ELSE 
      DO k=1, ke
        DO i = i_west_s+1, i_east_s-1
          ! radiation velocity = 0 for pp!
          bd_south(i,k, var_pp) = pp(i, j_south_s,   k)
        END DO
      END DO
    END IF

    DO k=1, ke
      DO i = i_west_s+1, i_east_s-1
        bd_south(i,k, var_Tp) = Tp(i, j_south_s,   k)                      &
          &         + cfl(1) * (Tp(i, j_south_s+1, k) - Tp(i, j_south_s, k) )
      END DO
    END DO

  END IF


  IF ( .NOT. l2dim ) THEN

    ! north-west corner
    IF ( (my_cart_neigh(2) == -1) .AND. (my_cart_neigh(1) == -1) ) THEN

      DO k=1, ke
        DO i = 1, i_west_u
          bd_north( i, k, var_u ) = u ( i_west_u+1, j_north_u-1, k)
        END DO
      END DO
      DO k=1, ke
        DO i = 1, i_west_v
          bd_north( i, k, var_v ) = v ( i_west_v+1, j_north_v-1, k)
        END DO
      END DO
      DO k=1, ke+1
        DO i = 1, i_west_s
          bd_north( i, k, var_w ) = w ( i_west_s+1, j_north_s-1, k)
        END DO
      END DO
      DO k=1, ke
        DO i = 1, i_west_s
          bd_north( i, k, var_pp) = pp( i_west_s+1, j_north_s-1, k)
          bd_north( i, k, var_Tp) = Tp( i_west_s+1, j_north_s-1, k)
        END DO
      END DO
    END IF

    ! north-east corner
    IF ( (my_cart_neigh(2) == -1) .AND. (my_cart_neigh(3) == -1)  ) THEN
      DO k=1, ke
        DO i = i_east_u, ie          
          bd_north( i, k, var_u ) = u ( i_east_u-1, j_north_u-1, k)
        END DO
      END DO
      DO k=1, ke
        DO i = i_east_v, ie
          bd_north( i, k, var_v ) = v ( i_east_v-1, j_north_v-1, k)
        END DO
      END DO
      DO k=1, ke+1
        DO i = i_east_s, ie
          bd_north( i, k, var_w ) = w ( i_east_s-1, j_north_s-1, k)
        END DO
      END DO
      DO k=1, ke
        DO i = i_east_s, ie
          bd_north( i, k, var_pp) = pp( i_east_s-1, j_north_s-1, k)
          bd_north( i, k, var_Tp) = Tp( i_east_s-1, j_north_s-1, k)
        END DO
      END DO

    END IF

    ! south-west corner
    IF ( (my_cart_neigh(4) == -1) .AND. (my_cart_neigh(1) == -1) ) THEN

      DO k=1, ke
        DO i = 1, i_west_u
          bd_south( i, k, var_u ) = u ( i_west_u+1, j_south_u+1, k)
        END DO
      END DO
      DO k=1, ke
        DO i = 1, i_west_v
          bd_south( i, k, var_v ) = v ( i_west_v+1, j_south_v+1, k)
        END DO
      END DO
      DO k=1, ke+1
        DO i = 1, i_west_s
          bd_south( i, k, var_w ) = w ( i_west_s+1, j_south_s+1, k)
        END DO
      END DO
      DO k=1, ke
        DO i = 1, i_west_s
          bd_south( i, k, var_pp) = pp( i_west_s+1, j_south_s+1, k)
          bd_south( i, k, var_Tp) = Tp( i_west_s+1, j_south_s+1, k)
        END DO
      END DO

    END IF

    ! south-east corner
    IF ( (my_cart_neigh(4) == -1) .AND. (my_cart_neigh(3) == -1) ) THEN

      DO k=1, ke
        DO i = i_east_u, ie
          bd_south( i, k, var_u ) = u ( i_east_u-1, j_south_u+1, k)
        END DO
      END DO
      DO k=1, ke
        DO i = i_east_v, ie
          bd_south( i, k, var_v ) = v ( i_east_v-1, j_south_v+1, k)
        END DO
      END DO
      DO k=1, ke+1
        DO i = i_east_s, ie
          bd_south( i, k, var_w ) = w ( i_east_s-1, j_south_s+1, k)
        END DO
      END DO
      DO k=1, ke
        DO i = i_east_s, ie
          bd_south( i, k, var_pp) = pp( i_east_s-1, j_south_s+1, k)
          bd_south( i, k, var_Tp) = Tp( i_east_s-1, j_south_s+1, k)
        END DO
      END DO

    END IF

  END IF

  DEALLOCATE( cfl )

END SUBROUTINE radiative_lbc

!=============================================================================
!=============================================================================

SUBROUTINE set_radiative_lbc( u, v, w, pp, Tp,                        &
  bd_west, bd_east, bd_north, bd_south )

  !--------------------------------------------------------------------------
  ! Description:
  !   Radiative lateral boundary condition for all dynamical variables.
  !   Here: set the values of u, v, w, pp, Tp at the boundary
  !   with the values bd_west, ... calculated in subroutine 'radiative_lbc'
  !
  ! Input variables:
  !   2-dim. boundary fields for the 5 variables:
  !      bd_west(:,:,5),  bd_east(:,:,5),  bd_north(:,:,5),  bd_south(:,:,5) 
  !
  ! Output variables:
  !   dynamical fields 
  !     u(:,:,:), v(:,:,:), w(:,:,:), pp(:,:,:), Tp(:,:,:)
  !   at the boundary points
  !
  !-------------------------- -----------------------------------------------

  REAL (kind=ireals), INTENT(in)  :: bd_west (je,ke1,5)
  REAL (kind=ireals), INTENT(in)  :: bd_east (je,ke1,5)
  REAL (kind=ireals), INTENT(in)  :: bd_north(ie,ke1,5)
  REAL (kind=ireals), INTENT(in)  :: bd_south(ie,ke1,5)

  REAL (kind=ireals), INTENT(inout) ::  u(ie,je,ke)
  REAL (kind=ireals), INTENT(inout) ::  v(ie,je,ke)
  REAL (kind=ireals), INTENT(inout) ::  w(ie,je,ke1)
  REAL (kind=ireals), INTENT(inout) :: pp(ie,je,ke)
  REAL (kind=ireals), INTENT(inout) :: Tp(ie,je,ke)

  INTEGER :: i, j, k

  ! Indizes for the position of u, v, s (=scalars pp, Tp, or w).
  INTEGER i_west_u
  INTEGER i_west_v
  INTEGER i_west_s
  INTEGER i_east_u
  INTEGER i_east_v
  INTEGER i_east_s
  INTEGER j_north_u
  INTEGER j_north_v
  INTEGER j_north_s
  INTEGER j_south_u
  INTEGER j_south_v
  INTEGER j_south_s

  ! these values must be the same as in subr. 'radiative_lbc':
  INTEGER :: var_u  = 1
  INTEGER :: var_v  = 2
  INTEGER :: var_w  = 3
  INTEGER :: var_Tp = 4
  INTEGER :: var_pp = 5


  WRITE(*,*) "subr. [set_radiative_lbc] ..."

  ! Indizes for the boundary position of u, v, s (=scalars pp, Tp, or w).
  ! They are the same as in subr. 'radiative_lbc':

  i_west_u = istartu-1
  i_west_v = istartv
  i_west_s = istart

  i_east_u = iendu+1
  i_east_v = iendv
  i_east_s = iend

  j_north_u = jendu
  j_north_v = jendv+1
  j_north_s = jend

  j_south_u = jstartu
  j_south_v = jstartv-1
  j_south_s = jstart


  ! western boundary
  IF ( my_cart_neigh(1) == -1 ) THEN

    DO k=1, ke
      DO j = j_south_u, j_north_u
        ! Anm.: eigentlich nur:  j = j_south_u+1, j_north_u-1,
        ! aber so hat man auch den Fall  l2dim=.TRUE.  dabei
        ! (analog folgende Schleifen)
        DO i = 1, i_west_u
          u (i, j, k) = bd_west(j, k, var_u ) 
        END DO
      END DO
    END DO

    IF ( .NOT. l2dim) THEN
      DO k=1, ke
        DO j = j_south_v+1, j_north_v-1
          DO i = 1, i_west_v
            v (i, j, k) = bd_west(j, k, var_v ) 
          END DO
        END DO
      END DO
    END IF

    DO k=1, ke1
      DO j = j_south_s, j_north_s
        DO i = 1, i_west_s
          w (i, j, k) = bd_west(j, k, var_w ) 
        END DO
      END DO
    END DO

    DO k=1, ke
      DO j = j_south_s, j_north_s
        DO i = 1, i_west_s
          pp(i, j, k) = bd_west(j, k, var_pp) 
          Tp(i, j, k) = bd_west(j, k, var_Tp) 
        END DO
      END DO
    END DO

  END IF

  ! eastern boundary
  IF ( my_cart_neigh(3) == -1 ) THEN

    DO k=1, ke
      DO j = j_south_u, j_north_u
        ! Anm.: eigentlich nur:  j = j_south_u+1, j_north_u-1,
        ! aber so hat man auch den Fall  l2dim=.TRUE.  dabei
        ! (analog folgende Schleifen)
        DO i = i_east_u, ie
          u (i, j, k) = bd_east(j, k, var_u ) 
        END DO
      END DO
    END DO

    IF ( .NOT. l2dim) THEN
      DO k=1, ke
        DO j = j_south_v+1, j_north_v-1
          DO i = i_east_v, ie 
            v (i, j, k) = bd_east(j, k, var_v ) 
          END DO
        END DO
      END DO
    END IF

    DO k=1, ke1
      DO j = j_south_s, j_north_s
        DO i = i_east_s, ie
          w (i, j, k) = bd_east(j, k, var_w ) 
        END DO
      END DO
    END DO

    DO k=1, ke
      DO j = j_south_s, j_north_s
        DO i = i_east_s, ie
          pp(i, j, k) = bd_east(j, k, var_pp ) 
          Tp(i, j, k) = bd_east(j, k, var_Tp ) 
        END DO
      END DO
    END DO

  END IF



  IF ( .NOT. l2dim ) THEN

    IF ( my_cart_neigh(2) == -1) THEN
      ! northern boundary (and NW-, NE-corners)

      DO k=1, ke
        DO j = j_north_u, je
          DO i = 1, ie
            u (i, j, k) = bd_north(i, k, var_u ) 
          END DO
        END DO
      END DO

      DO k=1, ke
        DO j = j_north_v, je
          DO i = 1, ie
            v (i, j, k) = bd_north(i, k, var_v ) 
          END DO
        END DO
      END DO

      DO k=1, ke1
        DO j = j_north_s, je
          DO i = 1, ie
            w (i, j, k) = bd_north(i, k, var_w ) 
          END DO
        END DO
      END DO

      DO k=1, ke
        DO j = j_north_s, je
          DO i = 1, ie
            pp(i, j, k) = bd_north(i, k, var_pp ) 
            Tp(i, j, k) = bd_north(i, k, var_Tp ) 
          END DO
        END DO
      END DO

    END IF


    ! southern boundary (and SW-, SE-corners)
    IF ( my_cart_neigh(4) == -1 ) THEN

      DO k=1, ke
        DO j = 1, j_south_u
          DO i = 1, ie
            u (i, j, k) = bd_south(i, k, var_u ) 
          END DO
        END DO
      END DO

      DO k=1, ke
        DO j = 1, j_south_v
          DO i = 1, ie
            v (i, j, k) = bd_south(i, k, var_v ) 
          END DO
        END DO
      END DO

      DO k=1, ke1
        DO j = 1, j_south_s
          DO i = 1, ie
            w (i, j, k) = bd_south(i, k, var_w ) 
          END DO
        END DO
      END DO

      DO k=1, ke
        DO j = 1, j_south_s
          DO i = 1, ie
            pp(i, j, k) = bd_south(i, k, var_pp ) 
            Tp(i, j, k) = bd_south(i, k, var_Tp ) 
          END DO
        END DO
      END DO

    END IF

  END IF

END SUBROUTINE set_radiative_lbc


!-----------------------------------------------------------------------------
!  End of the module fast_waves_sc
!-----------------------------------------------------------------------------

END MODULE fast_waves_sc

