#include "messy_main_ppd_bi.inc"

! **********************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL CONTRAIL
!
! Author : Volker Grewe, Christine Froemming, DLR-IPA, 2011
!
! References: see messy_contrail.f90
!
! **********************************************************************

! **********************************************************************
MODULE messy_contrail_si
! **********************************************************************

  ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      error_bi, info_bi

  ! SMCL
  USE messy_main_channel,       ONLY: t_chaobj_cpl
  USE messy_main_tools,         ONLY: PTR_3D_ARRAY, PTR_1D_ARRAY
  USE messy_contrail

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE
  SAVE

  ! GLOBAL PARAMETERS
  INTEGER, PARAMETER :: N_EMIS_MAX = 40

  ! GRID-POINT
  TYPE T_EMISSION_GP
     TYPE(t_chaobj_cpl) :: name
     REAL(DP), DIMENSION(:,:,:), POINTER :: emis => NULL()
     !
     ! OUTPUT
     REAL(DP), DIMENSION(:,:,:), POINTER :: concov  => NULL()  ! contrail cov.
     REAL(DP), DIMENSION(:,:,:), POINTER :: coniwc  => NULL()  ! contrail iwc
     REAL(DP), DIMENSION(:,:,:), POINTER :: coniwpc => NULL()  ! contrail iwp
     ! ...
  END TYPE T_EMISSION_GP
  TYPE(T_EMISSION_GP), DIMENSION(N_EMIS_MAX) :: XEMIS_GP
  INTEGER                                    :: N_EMIS_GP = 0

  ! op_sb_20140411+
  INTEGER, PARAMETER :: idp_cov = 1
  INTEGER, PARAMETER :: idp_ice = 2
  INTEGER, PARAMETER :: idp_max = 2
  CHARACTER(LEN=6), DIMENSION(2) :: pname = (/'concov','coniwc'/)
  ! op_sb_20140411-

  TYPE T_EMISSION_LG
     TYPE(t_chaobj_cpl) :: name
     REAL(DP), DIMENSION(:), POINTER :: emis => NULL()
     !
     ! OUTPUT
     REAL(DP), DIMENSION(:), POINTER :: coniwpc    => NULL() ! contrail iwp
     REAL(DP), DIMENSION(:), POINTER :: concov_sum => NULL() ! contrail coverage
     REAL(DP), DIMENSION(:), POINTER :: coniwc_sum => NULL() ! contrail iwc
     REAL(DP), DIMENSION(:), POINTER :: concov_now => NULL() ! contrail coverage
     REAL(DP), DIMENSION(:), POINTER :: coniwc_now => NULL() ! contrail iwc
     REAL(DP), DIMENSION(:), POINTER :: concov_m1  => NULL() ! contrail coverage
     REAL(DP), DIMENSION(:), POINTER :: coniwc_m1  => NULL() ! contrail iwc
     !
     REAL(DP), DIMENSION(:), POINTER :: te_spread => NULL() ! spread (tendency)
     REAL(DP), DIMENSION(:), POINTER :: te_sedi   => NULL() ! sedi (tendency)
     REAL(DP), DIMENSION(:), POINTER :: te_pot    => NULL() ! pot.cov (tendency)
     ! ...
     ! op_sb_20140411+
     TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER :: gp => NULL()
     ! op_sb_20140411-
     !
  END TYPE T_EMISSION_LG
  TYPE(T_EMISSION_LG), DIMENSION(N_EMIS_MAX) :: XEMIS_LG
  INTEGER                                    :: N_EMIS_LG = 0

  ! CPL-NAMELIST PARAMETERS
  LOGICAL            :: L_GP = .TRUE.
  TYPE(t_chaobj_cpl) :: C_GP_CLOUD_CRIT
  TYPE(t_chaobj_cpl) :: C_GP_CLOUD_COND
  TYPE(t_chaobj_cpl), DIMENSION(N_EMIS_MAX) :: C_GP_EMIS
  REAL(DP)           :: r_scal_gp = 1.0_dp
  !
  LOGICAL            :: L_LG = .FALSE.
  LOGICAL            :: L_LG_DIAG_TEND = .FALSE. ! diagnotic tendencies
  TYPE(t_chaobj_cpl) :: C_LG_CLOUD_CRIT
  TYPE(t_chaobj_cpl) :: C_LG_CLOUD_COND
  TYPE(t_chaobj_cpl), DIMENSION(N_EMIS_MAX) :: C_LG_EMIS
  REAL(DP)           :: r_scal_lg = 1.0_dp
  ! op_sb_20140411+
  ! calculate perturbed fields and transform to GP (LG contrail perturbation)
  LOGICAL            :: L_LG_calc_pert2GP = .FALSE.
  ! op_sb_20140411-

  ! POINTERS FOR COUPLED CHANNEL OBJECTS (BACKGROUND)
  REAL(DP), DIMENSION(:,:,:), POINTER :: gp_cloud_crit => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: gp_cloud_cond => NULL()
  !
  REAL(DP), DIMENSION(:),     POINTER :: lg_cloud_crit => NULL()
  REAL(DP), DIMENSION(:),     POINTER :: lg_cloud_cond => NULL()

  ! GOBAL PARAMETERS
  LOGICAL :: l_calc_lg = .FALSE.

  ! NEW CHANNEL OBJECTS: (DIAGNOSED FIELDS)
  ! GRID-POINT:
  ! ... POT. CONTRAIL COVERAGE
  ! 
  ! Thresholds for contrail cirrus
  REAL(DP), DIMENSION(:,:,:),   POINTER :: gp_b_cc => NULL()    ! pot. contr.-cir cov.
  REAL(DP), DIMENSION(:,:,:),   POINTER :: gp_potcov  => NULL() ! pot. contrail cov.
  REAL(DP), DIMENSION(:,:,:),   POINTER :: gp_qsm1 => NULL()    ! ?
!  REAL(DP), DIMENSION(:,:,:),   POINTER :: contr_qte => NULL()  ! 
!  REAL(DP), DIMENSION(:,:,:),   POINTER :: contr_aclc => NULL()  ! 
!  REAL(DP), DIMENSION(:,:,:),   POINTER :: contr_press => NULL()  ! 

  ! LAGRANGIAN:
  ! ... POT. CONTRAIL COVERAGE
  REAL(DP), DIMENSION(:),   POINTER :: lg_b_cc => NULL()    ! pot. contr.-cir cov.
  REAL(DP), DIMENSION(:),   POINTER :: lg_potcov  => NULL() ! pot. contrail cov.
  REAL(DP), DIMENSION(:),   POINTER :: lg_potcov_m1 => NULL() ! pot. contrail cov.
  REAL(DP), DIMENSION(:),   POINTER :: lg_qsm1 => NULL()    ! ?

! set pointer to GP cloud objects for transformation to lg, also if l_gp=false
  REAL(DP), DIMENSION(:,:,:),   POINTER :: cloud_crit => NULL()  ! 
  REAL(DP), DIMENSION(:,:,:),   POINTER :: cloud_cond => NULL()  ! 
! lg transformed from cloud 
  REAL(DP), DIMENSION(:),   POINTER :: cloud_crit_lg => NULL()
  REAL(DP), DIMENSION(:),   POINTER :: cloud_cond_lg => NULL()


  ! PUBLIC SUBROUTINES (called from messy_main_control_e5.f90)
  ! NOTE: in case you activate further entry points, make sure to call them
  !       in messy_main_control_e5.f90
  PUBLIC :: contrail_initialize    ! initialize submodel
  PUBLIC :: contrail_init_memory   ! request memory
  PUBLIC :: contrail_init_coupling ! initialize coupling to other SMs
  PUBLIC :: contrail_convec        ! op_lb_20160831
  PUBLIC :: contrail_physc         ! Calculates the physics
  PUBLIC :: contrail_global_end    ! Lagrangian

  ! PRIVATE SUBROTINES
  !PRIVATE :: contrail_read_nml_cpl

CONTAINS

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE contrail_initialize

    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,     ONLY: find_next_free_unit

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'contrail_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: i

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

    ! READ CTRL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find free I/O unit
       CALL contrail_read_nml_ctrl(status, iou)  ! read CTRL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CTRL namelist',substr)
    END IF
    ! BROADCAST CTRL namleist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(EI_H2O, p_io)
    CALL p_bcast(Q_fuel, p_io)
    CALL p_bcast(eta_ac,  p_io)
    CALL p_bcast(a_SAC,  p_io)
    CALL p_bcast(r_SAC,  p_io)

   ! READ CPL namelist
   IF (p_parallel_io) THEN                     ! read only on I/O-PE
      iou = find_next_free_unit(100,200)       ! find next free I/O unit
      CALL contrail_read_nml_cpl(status, iou)  ! read CPL-namelist
      ! terminate if error
      IF (status /= 0) CALL error_bi('Error in reading CPL namelist',substr)
   END IF
   ! BROADCAST CPL namleist entries from I/O-PE to ALL OTHER PEs
   CALL p_bcast(L_GP, p_io) 
   CALL p_bcast(r_scal_gp, p_io) 
   CALL p_bcast(C_GP_CLOUD_CRIT%cha, p_io) 
   CALL p_bcast(C_GP_CLOUD_CRIT%obj, p_io) 
   CALL p_bcast(C_GP_CLOUD_COND%cha, p_io) 
   CALL p_bcast(C_GP_CLOUD_COND%obj, p_io) 
   !
   CALL p_bcast(L_LG, p_io) 
   CALL p_bcast(r_scal_lg, p_io) 
   CALL p_bcast(L_LG_DIAG_TEND, p_io) 
   CALL p_bcast(l_calc_lg, p_io)             ! only set on I/O-PE
   CALL p_bcast(C_LG_CLOUD_CRIT%cha, p_io) 
   CALL p_bcast(C_LG_CLOUD_CRIT%obj, p_io) 
   CALL p_bcast(C_LG_CLOUD_COND%cha, p_io) 
   CALL p_bcast(C_LG_CLOUD_COND%obj, p_io) 
   CALL p_bcast(L_LG_calc_pert2GP,   p_io)   ! op_sb_20140411

   IF (L_GP) THEN
      IF (p_parallel_io) THEN
         N_EMIS_GP = 1
         DO i=1, N_EMIS_MAX
            IF (TRIM(C_GP_EMIS(i)%cha) == '') CYCLE
            IF (TRIM(C_GP_EMIS(i)%obj) == '') CYCLE
            XEMIS_GP(N_EMIS_GP)%name%cha = TRIM(C_GP_EMIS(i)%cha)
            XEMIS_GP(N_EMIS_GP)%name%obj = TRIM(C_GP_EMIS(i)%obj)
            ! NEXT EMISSION
            N_EMIS_GP = N_EMIS_GP + 1
         END DO
         N_EMIS_GP = N_EMIS_GP - 1
      END IF
      CALL p_bcast(N_EMIS_GP, p_io)
      !
      DO i=1, N_EMIS_GP
         CALL p_bcast(XEMIS_GP(i)%name%cha  , p_io)
         CALL p_bcast(XEMIS_GP(i)%name%obj  , p_io)
      END DO
   END IF

!!#D attila +
#ifdef ECHAM5
   IF (l_calc_lg) THEN
      IF (p_parallel_io) THEN
         N_EMIS_LG = 1
         DO i=1, N_EMIS_MAX
            IF (TRIM(C_LG_EMIS(i)%cha) == '') CYCLE
            IF (TRIM(C_LG_EMIS(i)%obj) == '') CYCLE
            XEMIS_LG(N_EMIS_LG)%name%cha = TRIM(C_LG_EMIS(i)%cha)
            XEMIS_LG(N_EMIS_LG)%name%obj = TRIM(C_LG_EMIS(i)%obj)
            ! NEXT EMISSION
            N_EMIS_LG = N_EMIS_LG + 1
         END DO
         N_EMIS_LG = N_EMIS_LG - 1
      END IF
      CALL p_bcast(N_EMIS_LG, p_io)
      !
      DO i=1, N_EMIS_LG
         CALL p_bcast(XEMIS_LG(i)%name%cha  , p_io)
         CALL p_bcast(XEMIS_LG(i)%name%obj  , p_io)
      END DO
   END IF
#endif
!!#D attila -

   CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

  END SUBROUTINE contrail_initialize
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE contrail_init_memory

    ! ------------------------------------------------------------------
    ! This subroutine is used to request memory for the submodel.
    ! The preferable method is to use "channel objects".
    ! Allocate your own memory, only if absolutely required.
    ! ------------------------------------------------------------------

    ! BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, LG_ATTILA
    USE messy_main_channel,          ONLY: new_channel, new_channel_object, &
                                           new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'contrail_init_memory'
    INTEGER                     :: status
    INTEGER                     :: i, j
    CHARACTER(LEN=3)            :: istr


    ! CHANNEL AND CHANNEL OBJECTS
    CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output

    gridpoint: IF (L_GP) THEN
       ! new channel
       CALL new_channel(status, modstr//'_gp', reprid=GP_3D_MID)
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//'_gp', 'b_cc', &
            p3=gp_b_cc)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_gp', 'b_cc', &
            'long_name', c='Potential contrail cirrus coverage')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_gp', 'b_cc', &
            'units', c='fraction')
       CALL channel_halt(substr, status)
       CALL info_bi('channel/object '//modstr//'_gp/'// &
            'b_cc'//' was created', substr)

       CALL new_channel_object(status, modstr//'_gp', 'potcov', &
            p3=gp_potcov)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_gp', 'potcov', &
            'long_name', c='potential contrail coverage')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_gp', 'potcov', &
            'units', c='fraction')
       CALL channel_halt(substr, status)
       CALL info_bi('channel/object '//modstr//'_gp/'// &
            'potcov'//' was created', substr)

       CALL new_channel_object(status, modstr//'_gp', 'qsm1', p3=gp_qsm1)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_gp', 'qsm1', 'long_name' &
            , c='qsm1')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_gp', 'qsm1', 'units', c='')
       CALL channel_halt(substr, status)
       CALL info_bi('channel/object '//modstr//'_gp/'// 'qsm1'//&
            &' was created', substr)

       DO i=1, N_EMIS_GP
          write(istr,'(i3.3)') i

          CALL new_channel_object(status, modstr//'_gp', 'concov_'//istr, &
               p3=XEMIS_GP(i)%concov)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_gp', 'concov_'//istr, &
               'long_name', c='contrail coverage')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_gp', 'concov_'//istr, &
               'units', c='fraction')
          CALL channel_halt(substr, status)
          CALL info_bi('channel/object '//modstr//'_gp/'// &
               'concov_'//istr//' was created', substr)

          CALL new_channel_object(status, modstr//'_gp', 'coniwc_'//istr, &
               p3=XEMIS_GP(i)%coniwc)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_gp', 'coniwc_'//istr, &
               'long_name', c='contrail ice water content')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_gp', 'coniwc_'//istr, &
               'units', c='kg/kg')
          CALL channel_halt(substr, status)
          CALL info_bi('channel/object '//modstr//'_gp/'// &
               'coniwc_'//istr//' was created', substr)

          CALL new_channel_object(status, modstr//'_gp', 'coniwpc_'//istr, &
               p3=XEMIS_GP(i)%coniwpc)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_gp', 'coniwpc_'//istr, &
               'long_name', c='contrail ice water path')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_gp', 'coniwpc_'//istr, &
               'units', c='g/m^2')
          CALL channel_halt(substr, status)
          CALL info_bi('channel/object '//modstr//'_gp/'// &
               'coniwpc_'//istr//' was created', substr)

       END DO

    END IF gridpoint

!!#D attila +
#ifdef ECHAM5
    Lagrangian: IF (l_calc_lg) THEN
       ! new channel
       CALL new_channel(status, modstr//'_lg', reprid=LG_ATTILA)
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//'_lg', 'b_cc', &
            p1=lg_b_cc)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'b_cc', &
            'long_name', c='Potential contrail cirrus coverage')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'b_cc', &
            'units', c='fraction')
       CALL channel_halt(substr, status)
       CALL info_bi('channel/object '//modstr//'_lg/'// &
            'b_cc'//' was created', substr)

       CALL new_channel_object(status, modstr//'_lg', 'potcov', &
            p1=lg_potcov)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'potcov', &
            'long_name', c='potential contrail coverage')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'potcov', &
            'units', c='fraction')
       CALL channel_halt(substr, status)
       CALL info_bi('channel/object '//modstr//'_lg/'// &
            'potcov'//' was created', substr)

       CALL new_channel_object(status, modstr//'_lg', 'potcov_m1', &
            p1=lg_potcov_m1, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'potcov_m1', &
            'long_name', c='potential contrail coverage at t-1')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'potcov_m1', &
            'units', c='fraction')
       CALL channel_halt(substr, status)
       CALL info_bi('channel/object '//modstr//'_lg/'// &
            'potcov_m1'//' was created', substr)

       CALL new_channel_object(status, modstr//'_lg', 'qsm1', p1=lg_qsm1)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'qsm1', 'long_name' &
            , c='qsm1')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'qsm1', 'units', c='')
       CALL channel_halt(substr, status)
       CALL info_bi('channel/object '//modstr//'_lg/'// 'qsm1'//&
            &' was created', substr)

! objects to save variables in contrail_physc for lg calculation in contrail_global_end
!       CALL new_channel_object(status, modstr//'_lg', 'contr_qte', &
!            p3=contr_qte, reprid=GP_3D_MID)
!       CALL channel_halt(substr, status)
!       CALL new_attribute(status, modstr//'_lg', 'contr_qte', &
!            'long_name', c='qte saved in contr_physc')
!       CALL channel_halt(substr, status)
!       CALL new_attribute(status, modstr//'_lg', 'contr_qte', &
!            'units', c='')
!       CALL channel_halt(substr, status)
!       CALL info_bi('channel/object '//modstr//'_lg/'// &
!            'contr_qte'//' was created', substr)
!
!       CALL new_channel_object(status, modstr//'_lg', 'contr_aclc', &
!            p3=contr_aclc, reprid=GP_3D_MID)
!       CALL channel_halt(substr, status)
!       CALL new_attribute(status, modstr//'_lg', 'contr_aclc', &
!            'long_name', c='aclc saved in contr_physc')
!       CALL channel_halt(substr, status)
!       CALL new_attribute(status, modstr//'_lg', 'contr_aclc', &
!            'units', c='')
!       CALL channel_halt(substr, status)
!       CALL info_bi('channel/object '//modstr//'_lg/'// &
!            'contr_aclc'//' was created', substr)

!       CALL new_channel_object(status, modstr//'_lg', 'cloud_crit', &
!            p3=cloud_crit, reprid=GP_3D_MID)
!       CALL channel_halt(substr, status)
!       CALL new_attribute(status, modstr//'_lg', 'cloud_crit', &
!            'long_name', c='rhc saved in contr_physc')
!       CALL channel_halt(substr, status)
!       CALL new_attribute(status, modstr//'_lg', 'cloud_crit', &
!            'units', c='')
!       CALL channel_halt(substr, status)
!       CALL info_bi('channel/object '//modstr//'_lg/'// &
!            'cloud_crit'//' was created', substr)

!       CALL new_channel_object(status, modstr//'_lg', 'cloud_cond', &
!            p3=cloud_cond, reprid=GP_3D_MID)
!       CALL channel_halt(substr, status)
!       CALL new_attribute(status, modstr//'_lg', 'cloud_cond', &
!            'long_name', c='cond saved in contr_physc')
!       CALL channel_halt(substr, status)
!       CALL new_attribute(status, modstr//'_lg', 'cloud_cond', &
!            'units', c='')
!       CALL channel_halt(substr, status)
!       CALL info_bi('channel/object '//modstr//'_lg/'// &
!            'cloud_cond'//' was created', substr)


!       CALL new_channel_object(status, modstr//'_lg', 'contr_press', &
!            p3=contr_press, reprid=GP_3D_MID)
!       CALL channel_halt(substr, status)
!       CALL new_attribute(status, modstr//'_lg', 'contr_press', &
!            'long_name', c='press saved in contr_physc')
!       CALL channel_halt(substr, status)
!       CALL new_attribute(status, modstr//'_lg', 'contr_press', &
!            'units', c='')
!       CALL channel_halt(substr, status)
!       CALL info_bi('channel/object '//modstr//'_lg/'// &
!            'contr_press'//' was created', substr)
! 
       CALL new_channel_object(status, modstr//'_lg', 'cloud_crit_lg', p1=cloud_crit_lg)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'cloud_crit_lg', 'long_name' &
            , c='cloud_crit_lg')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'cloud_crit_lg', 'units', c='')
       CALL channel_halt(substr, status)
       CALL info_bi('channel/object '//modstr//'_lg/'// 'cloud_crit_lg'//&
            &' was created', substr)


       CALL new_channel_object(status, modstr//'_lg', 'cloud_cond_lg', p1=cloud_cond_lg)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'cloud_cond_lg', 'long_name' &
            , c='cloud_cond_lg')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'cloud_cond_lg', 'units', c='')
       CALL channel_halt(substr, status)
       CALL info_bi('channel/object '//modstr//'_lg/'// 'cloud_cond_lg'//&
            &' was created', substr)


       DO i=1, N_EMIS_LG
          write(istr,'(i3.3)') i

          CALL new_channel_object(status, modstr//'_lg', 'concov_sum_'//istr, &
               p1=XEMIS_LG(i)%concov_sum)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'concov_sum_'//istr, &
               'long_name', c='contrail coverage (aged + new)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'concov_sum_'//istr, &
               'units', c='fraction')
          CALL channel_halt(substr, status)
          CALL info_bi('channel/object '//modstr//'_lg/'// &
               'concov_sum_'//istr//' was created', substr)

          CALL new_channel_object(status, modstr//'_lg', 'coniwc_sum_'//istr, &
               p1=XEMIS_LG(i)%coniwc_sum)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'coniwc_sum_'//istr, &
               'long_name', c='contrail ice water content  (aged + new)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'coniwc_sum_'//istr, &
               'units', c='kg/kg')
          CALL channel_halt(substr, status)
          CALL info_bi('channel/object '//modstr//'_lg/'// &
               'coniwc_sum_'//istr//' was created', substr)

          CALL new_channel_object(status, modstr//'_lg', 'concov_now_'//istr, &
               p1=XEMIS_LG(i)%concov_now)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'concov_now_'//istr, &
               'long_name', c='contrail coverage (actual emission)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'concov_now_'//istr, &
               'units', c='fraction')
          CALL channel_halt(substr, status)
          CALL info_bi('channel/object '//modstr//'_lg/'// &
               'concov_now_'//istr//' was created', substr)

          CALL new_channel_object(status, modstr//'_lg', 'coniwc_now_'//istr, &
               p1=XEMIS_LG(i)%coniwc_now)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'coniwc_now_'//istr, &
               'long_name', c='contrail ice water content (actual emission)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'coniwc_now_'//istr, &
               'units', c='kg/kg')
          CALL channel_halt(substr, status)
          CALL info_bi('channel/object '//modstr//'_lg/'// &
               'coniwc_now_'//istr//' was created', substr)

          CALL new_channel_object(status, modstr//'_lg', 'concov_m1_'//istr, &
               p1=XEMIS_LG(i)%concov_m1)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'concov_m1_'//istr, &
               'long_name', c='contrail coverage (at t-1)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'concov_m1_'//istr, &
               'units', c='fraction')
          CALL channel_halt(substr, status)
          CALL info_bi('channel/object '//modstr//'_lg/'// &
               'concov_m1_'//istr//' was created', substr)

          CALL new_channel_object(status, modstr//'_lg', 'coniwc_m1_'//istr, &
               p1=XEMIS_LG(i)%coniwc_m1)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'coniwc_m1_'//istr, &
               'long_name', c='contrail ice water content (at t-1)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'coniwc_m1_'//istr, &
               'units', c='kg/kg')
          CALL channel_halt(substr, status)
          CALL info_bi('channel/object '//modstr//'_lg/'// &
               'coniwc_m1_'//istr//' was created', substr)

          CALL new_channel_object(status, modstr//'_lg', 'coniwpc_'//istr, &
               p1=XEMIS_LG(i)%coniwpc)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'coniwpc_'//istr, &
               'long_name', c='contrail ice water path')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'coniwpc_'//istr, &
               'units', c='g/m^2')
          CALL channel_halt(substr, status)
          CALL info_bi('channel/object '//modstr//'_lg/'// &
               'coniwpc_'//istr//' was created', substr)

          IF (L_LG_DIAG_TEND) THEN

             CALL new_channel_object(status, modstr//'_lg', &
                  'te_spread_'//istr, &
                  p1=XEMIS_LG(i)%te_spread)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_lg', 'te_spread_'//istr, &
                  'long_name', c='contrail spreading tendency')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_lg', 'te_spread_'//istr, &
                  'units', c='1/s')
             CALL channel_halt(substr, status)
             CALL info_bi('channel/object '//modstr//'_lg/'// &
                  'te_spread_'//istr//' was created', substr)

             CALL new_channel_object(status, modstr//'_lg', &
                  'te_sedi_'//istr, &
                  p1=XEMIS_LG(i)%te_sedi)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_lg', 'te_sedi_'//istr, &
                  'long_name', c='tendency of contrail coverage by sedimentation')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_lg', 'te_sedi_'//istr, &
                  'units', c='kg/kg/s')
             CALL channel_halt(substr, status)
             CALL info_bi('channel/object '//modstr//'_lg/'// &
                  'te_sedi_'//istr//' was created', substr)

             CALL new_channel_object(status, modstr//'_lg', &
                  'te_pot_'//istr, &
                  p1=XEMIS_LG(i)%te_pot)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_lg', 'te_pot_'//istr, &
                  'long_name', c='tendency of potential contrail coverage')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_lg', 'te_pot_'//istr, &
                  'units', c='1/s')
             CALL channel_halt(substr, status)
             CALL info_bi('channel/object '//modstr//'_lg/'// &
                  'te_pot_'//istr//' was created', substr)

          END IF

          ! op_sb_20140411+
          IF (L_LG_calc_pert2GP) THEN
             ALLOCATE(XEMIS_LG(i)%gp(idp_max))
             DO j=1, idp_max
                CALL new_channel_object(status, modstr//'_lg', &
                     TRIM(pname(j))//'_'//istr,                &
                     p3 = XEMIS_LG(i)%gp(j)%ptr,             &
                     reprid = GP_3D_MID)
             END DO
          END IF
          ! op_sb_20140411-

       END DO

    END IF Lagrangian
#endif
!!#D attila -

    CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output

  END SUBROUTINE contrail_init_memory
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE contrail_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! MESSy
    USE messy_main_channel,          ONLY: get_channel_object

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'contrail_init_coupling'
    INTEGER :: i
    INTEGER :: status

    CALL start_message_bi(modstr,'COUPLING INITIALIZATION',substr)

    gridpoint: IF (L_GP) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) 'Grid-point calculation:'
          WRITE(*,*) 'Checking for ...'
          WRITE(*,*) '    channel: ',TRIM(C_GP_CLOUD_CRIT%cha)
          WRITE(*,*) '    object : ',TRIM(C_GP_CLOUD_CRIT%obj)
       END IF
       CALL get_channel_object(status, TRIM(C_GP_CLOUD_CRIT%cha) &
            , TRIM(C_GP_CLOUD_CRIT%obj)  &
            , p3=gp_cloud_crit)
       CALL channel_halt(substr, status)

       IF (p_parallel_io) THEN
          WRITE(*,*) 'Checking for ...'
          WRITE(*,*) '    channel: ',TRIM(C_GP_CLOUD_COND%cha)
          WRITE(*,*) '    object : ',TRIM(C_GP_CLOUD_COND%obj)
       END IF
       CALL get_channel_object(status, TRIM(C_GP_CLOUD_COND%cha) &
            , TRIM(C_GP_CLOUD_COND%obj)  &
            , p3=gp_cloud_cond)
       CALL channel_halt(substr, status)

       DO i=1, N_EMIS_GP
          IF (p_parallel_io) THEN
             WRITE(*,*) 'Checking for ...'
             WRITE(*,*) '    channel: ',TRIM(XEMIS_GP(i)%name%cha)
             WRITE(*,*) '    object : ',TRIM(XEMIS_GP(i)%name%obj)
          END IF
          CALL get_channel_object(status, TRIM(XEMIS_GP(i)%name%cha) &
               , TRIM(XEMIS_GP(i)%name%obj)  &
               , p3=XEMIS_GP(i)%emis)
          CALL channel_halt(substr, status)         
       END DO       
    END IF gridpoint

!!#D attila +
#ifdef ECHAM5
    Lagrangian: IF (l_calc_lg) THEN
!
       IF (p_parallel_io) THEN
          WRITE(*,*) 'Lagrangian calculation:'
          WRITE(*,*) 'Checking for ...'
          WRITE(*,*) '    channel: cloud '
          WRITE(*,*) '    object : rhc '
       END IF
       CALL get_channel_object(status,'cloud' &
            , 'rhc'   &
            , p3=cloud_crit)
       CALL channel_halt(substr, status)

       IF (p_parallel_io) THEN
          WRITE(*,*) 'Lagrangian calculation:'
          WRITE(*,*) 'Checking for ...'
          WRITE(*,*) '    channel: cloud '
          WRITE(*,*) '    object : condensation '
       END IF
       CALL get_channel_object(status,'cloud' &
            , 'condensation'   &
            , p3=cloud_cond)
       CALL channel_halt(substr, status)


       IF (p_parallel_io) THEN
          WRITE(*,*) 'Lagrangian calculation:'
          WRITE(*,*) 'Checking for ...'
          WRITE(*,*) '    channel: ',TRIM(C_LG_CLOUD_CRIT%cha)
          WRITE(*,*) '    object : ',TRIM(C_LG_CLOUD_CRIT%obj)
       END IF
       CALL get_channel_object(status, TRIM(C_LG_CLOUD_CRIT%cha) &
            , TRIM(C_LG_CLOUD_CRIT%obj)  &
            , p1=lg_cloud_crit)
       CALL channel_halt(substr, status)

       IF (p_parallel_io) THEN
          WRITE(*,*) 'Checking for ...'
          WRITE(*,*) '    channel: ',TRIM(C_LG_CLOUD_COND%cha)
          WRITE(*,*) '    object : ',TRIM(C_LG_CLOUD_COND%obj)
       END IF
       CALL get_channel_object(status, TRIM(C_LG_CLOUD_COND%cha) &
            , TRIM(C_LG_CLOUD_COND%obj)  &
            , p1=lg_cloud_cond)
       CALL channel_halt(substr, status)

       DO i=1, N_EMIS_LG
          IF (p_parallel_io) THEN
             WRITE(*,*) 'Checking for ...'
             WRITE(*,*) '    channel: ',TRIM(XEMIS_LG(i)%name%cha)
             WRITE(*,*) '    object : ',TRIM(XEMIS_LG(i)%name%obj)
          END IF
          CALL get_channel_object(status, TRIM(XEMIS_LG(i)%name%cha) &
               , TRIM(XEMIS_LG(i)%name%obj)  &
               , p1=XEMIS_LG(i)%emis)
          CALL channel_halt(substr, status)         
       END DO
    END IF Lagrangian
#endif
!!#D attila -

    CALL end_message_bi(modstr,'COUPLING INITIALIZATION',substr)

  END SUBROUTINE contrail_init_coupling
  ! ====================================================================

  ! ====================================================================
  ! op_lb_20160831+
  SUBROUTINE contrail_convec

    IMPLICIT NONE

    CALL contrail_physc

  END SUBROUTINE contrail_convec
  ! op_lb_20160831-
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE contrail_physc

    ! ------------------------------------------------------------------
    ! This subroutine is called within the time loop.
    ! It constitutes the main entry point for additional processes 
    ! or diagnostics.
    ! Here, only the current vector of the grid-point-fields is
    ! accessible.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_timer,           ONLY: time_step_len
    USE messy_main_grid_def_mem_bi, ONLY: nproma, kproma, jrow, nlev &
         , ngpblks
    USE messy_main_grid_def_bi,     ONLY: gboxarea_2d
    USE messy_main_data_bi,         ONLY:&
           press_3d    & ! mid-level pressures [Pa]
         , pressi_3d   &
         , tm1         &
         , qm1         &
         , qte_3d      &
         , aclc
    USE messy_main_constants_mem, ONLY: g

    IMPLICIT NONE

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr = 'contrail_physc'
    INTEGER                     :: jk, i
    REAL(dp)                    :: zdp(nproma,nlev)
    REAL(dp)                    :: zconpn(nproma,nlev)
    REAL(DP), DIMENSION(:,:), POINTER :: gboxarea  => NULL()

! save variables for lg_calculation

!    contr_qte(_RI_XYZ__(1:kproma,jrow,:))=qte_3d(_RI_XYZ__(1:kproma,jrow,:))
!    contr_aclc(_RI_XYZ__(1:kproma,jrow,:))=aclc(_RI_XYZ__(1:kproma,jrow,:))
!    contr_press(_RI_XYZ__(1:kproma,jrow,:))=press_3d(_RI_XYZ__(1:kproma,jrow,:))
! TO DO QQQ
!    cloud_crit(_RI_XYZ__(1:kproma,jrow,:))=gp_cloud_crit(_RI_XYZ__(1:kproma,jrow,:))
!    cloud_cond(_RI_XYZ__(1:kproma,jrow,:))=gp_cloud_cond(_RI_XYZ__(1:kproma,jrow,:))

    IF (.NOT. L_GP) RETURN

!qqq NOTE: use correct start values (x = x_m1 + x_tte * Dt ???)
!          for temp, q

    CALL contrail_pot_cov( &
         tm1(_RI_XYZ__(1:kproma,jrow,:)),           & ! IN:  temp
         press_3d(_RI_XYZ__(1:kproma,jrow,:)),      & ! IN:  press
         qm1(_RI_XYZ__(1:kproma,jrow,:)),           & ! IN:  water vapour
         qte_3d(_RI_XYZ__(1:kproma,jrow,:)),        & ! IN:  water vapour tend. from cloud
         aclc(_RI_XYZ__(1:kproma,jrow,:)),          & ! IN:  coverage
         gp_cloud_crit(_RI_XYZ__(1:kproma,jrow,:)), & ! IN:  critical hum. for nat. clouds
         time_step_len,                   & ! IN 
         gp_b_cc(_RI_XYZ__(1:kproma,jrow,:)),       & ! OUT: pot. contrail cirrus coverage
         gp_potcov(_RI_XYZ__(1:kproma,jrow,:)),     & ! OUT: potential contrail coverage
         gp_qsm1(_RI_XYZ__(1:kproma,jrow,:))  ,     & ! 
         zconpn(1:kproma,:) )

    IF (N_EMIS_GP > 0) THEN

       ALLOCATE(gboxarea(nproma, nlev))

       DO i=1,nlev
          gboxarea(1:kproma,i) = gboxarea_2d(1:kproma,jrow)
       END DO

       zdp(1:kproma,1:nlev) = pressi_3d(_RI_XYZ__(1:kproma,jrow,2:nlev+1)) -  &
            pressi_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))
    END IF

    DO i=1, N_EMIS_GP

       CALL contrail_calc( &
            qm1(_RI_XYZ__(1:kproma,jrow,:)),                & ! IN: water vapour
            aclc(_RI_XYZ__(1:kproma,jrow,:)),               & ! IN: coverage
            gp_cloud_cond(_RI_XYZ__(1:kproma,jrow,:)),      & ! IN: cond. rate from cloud
            gp_potcov(_RI_XYZ__(1:kproma,jrow,:)),          & ! IN: pot. contrail coverage
            zconpn(1:kproma,:),                   & ! IN:
            XEMIS_GP(i)%emis(_RI_XYZ__(1:kproma,jrow,:)),   & ! IN: emission
            r_scal_gp,                            & ! IN: scaling factor
            gboxarea(1:kproma,:),                           & ! IN: boxarea
            XEMIS_GP(i)%concov(_RI_XYZ__(1:kproma,jrow,:)), & ! OUT: contrail coverage
            XEMIS_GP(i)%coniwc(_RI_XYZ__(1:kproma,jrow,:))  & ! OUT: contr. ice water cont
            )

       ! factor 1000 for conversion from kg to g
       ! iwpc is ice water path (g/m^2) 
       XEMIS_GP(i)%coniwpc(_RI_XYZ__(1:kproma,jrow,:)) = &
            (XEMIS_GP(i)%coniwc(_RI_XYZ__(1:kproma,jrow,:)) &
            * 1000._dp * zdp(1:kproma,:)) / g

    END DO

    IF (ASSOCIATED(gboxarea)) THEN
       DEALLOCATE(gboxarea)
       NULLIFY(gboxarea)
    END IF

  END SUBROUTINE contrail_physc
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE contrail_global_end

!!#D attila +
#ifdef ECHAM5
    USE messy_main_timer,         ONLY: time_step_len
    USE messy_main_tracer_mem_bi, ONLY: NCELL
    USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev, npromz, ngpblks
    USE messy_main_grid_def_bi,     ONLY: sqcst_2d, gboxarea_2d
    USE messy_main_data_bi,         ONLY:  &
           press_3d    & ! mid-level pressures [Pa]
         , pressi_3d   &
         , etadot_3d   &
         , tm1         &
         , qm1         &
         , xim1_3d     &     ! op_sb_20140414
         , qte_3d      &
         , aclc        &
         , u_scb, v_scb
    USE messy_attila_tools_e5,    ONLY: gp2lg_e5,         &
                                        lg2gp_e5, LG2GP_AVE   ! op_sb_20140414
    USE messy_main_constants_mem, ONLY: g, rd

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! LOCAL
    INTEGER  :: kproma, jg, i, jk
    REAL(DP), DIMENSION(:), POINTER :: lg_zconpn => NULL()
    REAL(DP), DIMENSION(:), POINTER :: lg_tm1   => NULL()
    REAL(DP), DIMENSION(:), POINTER :: lg_press => NULL()
    REAL(DP), DIMENSION(:), POINTER :: lg_qm1   => NULL()
    REAL(DP), DIMENSION(:), POINTER :: lg_qte   => NULL()
    REAL(DP), DIMENSION(:), POINTER :: lg_aclc  => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: gp_zdp  => NULL()
    REAL(DP), DIMENSION(:),     POINTER :: lg_zdp  => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: gboxarea  => NULL()
    REAL(DP), DIMENSION(:),     POINTER :: lg_gboxarea  => NULL()
    REAL(DP), DIMENSION(:),     POINTER :: lg_rho => NULL()
    REAL(DP), DIMENSION(:),     POINTER :: lg_etadot => NULL()
    REAL(DP), DIMENSION(:),     POINTER :: lg_dudz => NULL()
    REAL(DP), DIMENSION(:),     POINTER :: lg_dvdz => NULL()
    !
    REAL(DP), DIMENSION(NCELL) :: te_spread
    REAL(DP), DIMENSION(NCELL) :: te_sedi
    REAL(DP), DIMENSION(NCELL) :: te_pot
    !
    REAL(DP), DIMENSION(:,:,:), POINTER :: rho => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: dudz => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: dvdz => NULL()
    
    IF (.NOT. l_calc_lg) RETURN

    ALLOCATE(lg_zconpn(NCELL))
    !
    ALLOCATE(lg_tm1(NCELL))
    ALLOCATE(lg_press(NCELL))
    ALLOCATE(lg_qm1(NCELL))
    ALLOCATE(lg_qte(NCELL))
    ALLOCATE(lg_aclc(NCELL))

    CALL gp2lg_e5(tm1,lg_tm1)
    CALL gp2lg_e5(qm1,lg_qm1)
    CALL gp2lg_e5(qte_3d,lg_qte)
    CALL gp2lg_e5(aclc,lg_aclc)
    CALL gp2lg_e5(press_3d,lg_press)
!    CALL gp2lg_e5(contr_press,lg_press)
!    CALL gp2lg_e5(contr_qte,lg_qte)
!    CALL gp2lg_e5(contr_aclc,lg_aclc)
    CALL gp2lg_e5(cloud_crit,cloud_crit_lg)
    CALL gp2lg_e5(cloud_cond,cloud_cond_lg)

    CALL contrail_pot_cov( &
         lg_tm1(:),        & ! IN:  temp
         lg_press(:),      & ! IN:  press
         lg_qm1(:),        & ! IN:  water vapour
         lg_qte(:),        & ! IN:  water vapour tendency from cloud
         lg_aclc(:),       & ! IN:  coverage
         lg_cloud_crit(:), & ! IN:  critical humidity for nat. clouds
         time_step_len,    & ! IN 
         lg_b_cc(:),       & ! OUT: potential contrail cirrus coverage
         lg_potcov(:),     & ! OUT: potential contrail coverage
         lg_qsm1(:),       &
         lg_zconpn(:)      )

    emissions: IF (N_EMIS_LG > 0) THEN

       ALLOCATE(gp_zdp(nproma, nlev, ngpblks))
       ALLOCATE(gboxarea(nproma, nlev, ngpblks))
       ALLOCATE(rho(nproma, nlev, ngpblks))
       ALLOCATE(dudz(nproma, nlev, ngpblks))
       ALLOCATE(dvdz(nproma, nlev, ngpblks))
       ALLOCATE(lg_zdp(NCELL))
       ALLOCATE(lg_gboxarea(NCELL))
       ALLOCATE(lg_etadot(NCELL))
       ALLOCATE(lg_rho(NCELL))
       ALLOCATE(lg_dudz(NCELL))
       ALLOCATE(lg_dvdz(NCELL))

       gp_zdp(:,1:nlev,:) = pressi_3d(:,2:nlev+1,:) - pressi_3d(:,1:nlev,:)
       call gp2lg_e5(gp_zdp,lg_zdp)

       do jk=1,nlev 
          gboxarea(:,jk,:) = gboxarea_2d(:,:) 
       end do
       call gp2lg_e5(gboxarea,lg_gboxarea)

       CALL gp2lg_e5(etadot_3d,lg_etadot)

       do jg = 1, ngpblks
          if (jg == ngpblks) then
             kproma = npromz
          else
             kproma = nproma
          endif
          rho(1:kproma,:,jg) = press_3d(1:kproma,:,jg)/(rd*tm1(1:kproma,:,jg))

          CALL contrail_uv_grad(u_scb(1:kproma,:,jg) &
               , v_scb(1:kproma,:,jg)    &
               , sqcst_2d(1:kproma,jg)   &
               , rho(1:kproma,:,jg)      &
               , press_3d(1:kproma,:,jg) &
               , dudz(1:kproma,:,jg)     &
               , dvdz(1:kproma,:,jg)     &
               )
       enddo
       call gp2lg_e5(rho,lg_rho)
       call gp2lg_e5(dudz,lg_dudz)
       call gp2lg_e5(dvdz,lg_dvdz)

    END IF emissions

    DO i=1, N_EMIS_LG

       CALL contrail_calc_dev( &
            lg_qm1(:),         & ! IN: water vapour
            lg_aclc(:),        & ! IN: coverage
            lg_cloud_cond(:),  & ! IN: condensation rate from cloud
            lg_potcov(:),      & ! IN: potential contrail coverage
            lg_gboxarea(:),    & ! IN:
            lg_zconpn(:),      & ! IN:
            time_step_len,            & ! IN:
            XEMIS_LG(i)%emis(:),      & ! IN: emission
            r_scal_lg,                &
            lg_potcov_m1(:),          & ! IN:
            XEMIS_LG(i)%concov_m1(:), &
            XEMIS_LG(i)%coniwc_m1(:), &
            lg_rho(:), lg_etadot(:),  &
            lg_dudz(:), lg_dvdz(:),   & ! IN     
            XEMIS_LG(i)%concov_sum(:), & ! OUT: contrail coverage
            XEMIS_LG(i)%coniwc_sum(:), & ! OUT: contrail ice water content
            XEMIS_LG(i)%concov_now(:), & ! OUT: contrail coverage
            XEMIS_LG(i)%coniwc_now(:), & ! OUT: contrail ice water content
            te_spread(:),  &
            te_sedi(:),    &
            te_pot(:)      &
            )

   ! SAVE FOR NEXT TIME STEP
          XEMIS_LG(i)%concov_m1(:) = XEMIS_LG(i)%concov_sum(:)
          XEMIS_LG(i)%coniwc_m1(:) = XEMIS_LG(i)%coniwc_sum(:)


       IF (L_LG_DIAG_TEND) THEN
          XEMIS_LG(i)%te_spread(:) = te_spread(:)
          XEMIS_LG(i)%te_sedi(:)   = te_sedi(:)
          XEMIS_LG(i)%te_pot(:)    = te_pot(:)
       END IF

       ! op_sb_20140414+
       ! make contrail cover and ice available on ECHAM grid
       IF (L_LG_calc_pert2GP) THEN
          CALL lg2gp_e5(XEMIS_LG(i)%concov_sum, XEMIS_LG(i)%gp(idp_cov)%ptr &
               , LG2GP_AVE)   ! contrail cover

          CALL lg2gp_e5(XEMIS_LG(i)%coniwc_sum, XEMIS_LG(i)%gp(idp_ice)%ptr &
               , LG2GP_AVE)   ! contrail icemixr

          DO jg = 1, ngpblks
             IF (jg == ngpblks) then
                kproma = npromz
             ELSE
                kproma = nproma
             END IF

          END DO
       END IF
       ! op_sb_20140414- 

    END DO

    ! SAVE FOR NEXT TIME STEP
    lg_potcov_m1(:) = lg_potcov(:)

    ! CLEAN UP MEMOPRY
    DEALLOCATE(lg_zconpn) ; NULLIFY(lg_zconpn)
    DEALLOCATE(lg_tm1)    ; NULLIFY(lg_tm1)
    DEALLOCATE(lg_press)  ; NULLIFY(lg_press)
    DEALLOCATE(lg_qm1)    ; NULLIFY(lg_qm1)
    DEALLOCATE(lg_qte)    ; NULLIFY(lg_qte)
    DEALLOCATE(lg_aclc)   ; NULLIFY(lg_aclc)

    IF (ASSOCIATED(gp_zdp)) THEN
       DEALLOCATE(gp_zdp)
       NULLIFY(gp_zdp)
    END IF

    IF (ASSOCIATED(gboxarea)) THEN
       DEALLOCATE(gboxarea)
       NULLIFY(gboxarea)
    END IF

    IF (ASSOCIATED(rho))  THEN
       DEALLOCATE(rho)
       NULLIFY(rho)
    END IF

    IF (ASSOCIATED(dudz)) THEN
       DEALLOCATE(dudz)
       NULLIFY(dudz)
    END IF

    IF (ASSOCIATED(dvdz)) THEN
       DEALLOCATE(dvdz)
       NULLIFY(dvdz)
    END IF

    IF (ASSOCIATED(lg_gboxarea)) THEN
       DEALLOCATE(lg_gboxarea)
       NULLIFY(lg_gboxarea)
    END IF

    IF (ASSOCIATED(lg_zdp)) THEN
       DEALLOCATE(lg_zdp)
       NULLIFY(lg_zdp)
    END IF

    IF (ASSOCIATED(lg_rho)) THEN
       DEALLOCATE(lg_rho)
       NULLIFY(lg_rho)
    END IF

    IF (ASSOCIATED(lg_dudz)) THEN
       DEALLOCATE(lg_dudz)
       NULLIFY(lg_dudz)
    END IF

    IF (ASSOCIATED(lg_dvdz)) THEN
       DEALLOCATE(lg_dvdz)
       NULLIFY(lg_dvdz)
    END IF

    IF (ASSOCIATED(lg_etadot)) THEN
       DEALLOCATE(lg_etadot)
       NULLIFY(lg_etadot)
    END IF

#endif
!!#D attila -

  END SUBROUTINE contrail_global_end
  ! ====================================================================

! ======================================================================
SUBROUTINE contrail_read_nml_cpl(status, iou)

  ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_main_tracer_mem_bi, ONLY: NGCELL

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ L_GP, C_GP_CLOUD_CRIT, C_GP_CLOUD_COND &
                 , L_LG, C_LG_CLOUD_CRIT, C_LG_CLOUD_COND &
                 , L_LG_DIAG_TEND, r_scal_gp, r_scal_lg   &
                 , L_LG_calc_pert2GP                      &  ! op_sb_20140604
                 , C_LG_EMIS, C_GP_EMIS                      ! op_cf_20140611

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='contrail_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
#ifdef ECHAM5
    IF ((L_LG) .AND. (NGCELL > 0)) THEN
       l_calc_lg = .TRUE.
!!#D attila +
       WRITE(*,*) 'Lagrangian: ON'
       IF (L_LG_calc_pert2GP) THEN
          WRITE(*,*) '... conversion of LG perturbation fields to GP (for CLOUDOPT) '
       ELSE
          WRITE(*,*) '... no conversion of LG perturbation fields to GP (for CLOUDOPT) '
       ENDIF
!!#D attila -
    ELSE
       IF (L_LG) THEN
!!#D attila +
         WRITE(*,*) 'L_LG = T in namelist'
         WRITE(*,*) 'However no Lagrangian scheme activated ...'
         WRITE(*,*) ' ... setting l_calc_lg = F'
!!#D attila -
       ENDIF
       l_calc_lg = .FALSE.
!!#D attila +
       WRITE(*,*) 'Lagrangian: OFF'
!!#D attila -
    ENDIF
#endif

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

END SUBROUTINE contrail_read_nml_cpl
! ======================================================================

! **********************************************************************
END MODULE messy_contrail_si
! **********************************************************************
