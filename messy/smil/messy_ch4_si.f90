! **********************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL CH4 
!
! Author : Patrick Joeckel, DLR-IPA, May  2012
!
! References: see messy_ch4.f90
!
! **********************************************************************
#include "messy_main_ppd_bi.inc"

! **********************************************************************
MODULE messy_ch4_si
! **********************************************************************

  ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      error_bi, warning_bi

  USE messy_main_channel,       ONLY: t_chaobj_cpl
  ! ISO++
  USE messy_main_constants_mem, ONLY: MC, MH, MO, M13C, M12C, MD
  ! ISO--
#ifdef MESSYTENDENCY
 !tendency budget
 USE messy_main_tendency_bi,    ONLY: mtend_get_handle,       &
                                      mtend_get_start_l,      &
                                      mtend_get_start_g,      &
                                      mtend_add_l,            &
                                      mtend_add_g,            &
                                      mtend_register,         &    
                                      mtend_id_tracer, mtend_id_q, mtend_id_t
#endif
  USE messy_main_timer_event,  ONLY: time_event, io_time_event

  ! SMCL
  USE messy_ch4

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE
  SAVE

  ! GLOBAL PARAMETERS
  INTEGER :: idt_gp_CH4 = 0    ! tracer ID (GP)
#ifdef ECHAM5
  INTEGER :: idt_lg_CH4 = 0    ! tracer ID (LG)
  INTEGER :: idt_lg_H2O = 0    ! tracer ID (LG)
#endif
#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  ! ISO++
  ! Array of used tracer ids
  INTEGER, DIMENSION(:), POINTER :: idt_gp_iso => NULL()
  INTEGER, DIMENSION(:), POINTER :: idt_lg_iso => NULL()

  REAL(DP),         DIMENSION(4), PARAMETER :: iso_molarmass = &
       (/M12C + 4.0_dp*MH, &               ! 12CH4
       M13C + 4.0_dp*MH, &                 ! 13CH4
       MC + 4.0_dp*MH, &                   ! CH4D0
       MC + 3.0_dp*MH + MD/)               ! CH3D1  

  REAL(DP), PARAMETER            :: hdo_molarmass = (MH + MD + MO)

  ! Array of used isotopologue names
  INTEGER, DIMENSION(:), POINTER :: iso_id_gp => NULL()
  INTEGER, DIMENSION(:), POINTER :: iso_id_lg => NULL()

  ! HDO++
  ! GP:
  INTEGER       :: i_gp_CH4_D1     = 0       ! index of CH3D tracer (GP)
  INTEGER       :: idt_gp_HDO      = 0       ! id of HDO tracer
  LOGICAL       :: HDO_ex_gp       = .FALSE. ! flag if HDO tracer exists
#ifdef ECHAM5
  ! LG:
  INTEGER       :: i_lg_CH4_D1     = 0       ! index of CH3D tracer (LG)
  INTEGER       :: idt_lg_HDO      = 0       ! id of HDO tracer
  LOGICAL       :: HDO_ex_lg       = .FALSE. ! flag if HDO tracer exists
#endif
  ! HDO--
  ! ISO--

  ! CPL-NAMELIST PARAMETERS
  INTEGER :: i_H2O_feedback = 0 ! no feedback to hydrological cycle
  !                             ! (1: feedback from GP; 2: feedback from LG)
  !
  LOGICAL :: L_GP = .TRUE.
  TYPE(t_chaobj_cpl) :: c_gp_OH
  TYPE(t_chaobj_cpl) :: c_gp_O1D
  TYPE(t_chaobj_cpl) :: c_gp_Cl
  TYPE(t_chaobj_cpl) :: c_gp_jCH4
  INTEGER, DIMENSION(2) :: i_gp_nclass_emis_age = (/0, 0/) ! op_pj_20140628
  !
  LOGICAL :: L_LG = .FALSE.
  TYPE(t_chaobj_cpl) :: c_lg_OH
  TYPE(t_chaobj_cpl) :: c_lg_O1D
  TYPE(t_chaobj_cpl) :: c_lg_Cl
  TYPE(t_chaobj_cpl) :: c_lg_jCH4
  ! CLASS++
  INTEGER, DIMENSION(2) :: i_lg_nclass_emis_age = (/0, 0/)
  LOGICAL :: l_gp_adj_tend = .TRUE.
  LOGICAL :: l_lg_adj_tend = .TRUE.
  INTEGER :: i_gp_ageing   = 1      ! ageing method (class transfer)
  INTEGER :: i_lg_ageing   = 1      ! ageing method (class transfer)
  ! CLASS--

  ! CLASS++
  LOGICAL :: l_gp_class = .FALSE.
  INTEGER,  DIMENSION(:,:),   POINTER :: idt_gp_class => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: sadj_gp => NULL()
  LOGICAL  :: l_gp_mv_agecl_now = .FALSE. ! move age classes now ?
  REAL(DP), DIMENSION(:), POINTER :: f_gp => NULL()
  !
  LOGICAL :: l_lg_class = .FALSE.
  INTEGER,  DIMENSION(:,:), POINTER :: idt_lg_class => NULL()
#ifdef ECHAM5
  REAL(DP), DIMENSION(:),   POINTER :: sadj_lg => NULL()
#endif
  LOGICAL  :: l_lg_mv_agecl_now = .FALSE. ! move age classes now ?
  REAL(DP), DIMENSION(:), POINTER :: f_lg => NULL()
  !
  INTEGER, PARAMETER :: NAGECLASSMAX = 20
  REAL(DP), DIMENSION(NAGECLASSMAX) :: r_gp_age_cll  
  REAL(DP), DIMENSION(NAGECLASSMAX) :: r_lg_age_cll  
  !
  ! TIMER
  ! only required for ageing method 0 (not recommended, see below!)
  TYPE(io_time_event) :: TIMER_MONTHLY = io_time_event(1, 'months','last',0)
  TYPE(time_event)    :: XTIMER_MONTHLY 
  ! CLASS--

  ! ISO++
  LOGICAL :: l_gp_iso_C = .FALSE. ! Isotopologues concerning C included
  LOGICAL :: l_gp_iso_H = .FALSE. ! Isotopologues concerning H included
  LOGICAL :: l_lg_iso_C = .FALSE. ! Isotopologues concerning C included
  LOGICAL :: l_lg_iso_H = .FALSE. ! Isotopologues concerning H included

  REAL(DP),  DIMENSION(:,:,:), POINTER :: sadj_iso_C_gp => NULL()
  REAL(DP),  DIMENSION(:,:,:), POINTER :: sadj_iso_H_gp => NULL()
#ifdef ECHAM5
  REAL(DP),  DIMENSION(:),     POINTER :: sadj_iso_C_lg => NULL()
  REAL(DP),  DIMENSION(:),     POINTER :: sadj_iso_H_lg => NULL()
#endif
  ! ISO--
  
  ! POINTERS
  REAL(dp), DIMENSION(:,:,:), POINTER :: ptr_gp_OH   => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: ptr_gp_O1D  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: ptr_gp_Cl   => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: ptr_gp_jCH4 => NULL()
  !
#ifdef ECHAM5
  REAL(dp), DIMENSION(:), POINTER :: ptr_lg_OH   => NULL()
  REAL(dp), DIMENSION(:), POINTER :: ptr_lg_O1D  => NULL()
  REAL(dp), DIMENSION(:), POINTER :: ptr_lg_Cl   => NULL()
  REAL(dp), DIMENSION(:), POINTER :: ptr_lg_jCH4 => NULL()
  !
  REAL(dp), DIMENSION(:), POINTER :: press_3d_lg => NULL()
#endif
  ! LOGICAL if empirical formula of roland eichinger should be used
  LOGICAL :: l_ef_re = .FALSE. 

  ! PUBLIC SUBROUTINES (called from messy_main_control_e5.f90)
  PUBLIC :: ch4_initialize    ! initialize submodel
  PUBLIC :: ch4_new_tracer    ! define new tracers
  PUBLIC :: ch4_init_memory   ! define channel objects
  PUBLIC :: ch4_init_coupling ! set pointers for coupling to BM and other SMs
  PUBLIC :: ch4_global_start
  PUBLIC :: ch4_vdiff         ! 
  PUBLIC :: ch4_physc         ! entry point in time loop (current vector)
  PUBLIC :: ch4_global_end    ! entry point in time loop for LG calculation
  PUBLIC :: ch4_free_memory   ! deallocate memory

  ! PRIVATE SUBROUTINES
  !PRIVATE :: ch4_read_nml_cpl
  !
  !CLASS++
  !PRIVATE :: class_integrate_gp  ! (called from ch4_physc)
  !PRIVATE :: class_age_move_gp   ! (called from class_integrate_gp)
  !PRIVATE :: class_adj_tend_gp   ! (called from class_integrate_gp)
  !
  !PRIVATE :: class_integrate_lg ! (called from ch4_global_end)
  !PRIVATE :: class_age_move_lg  ! (called from class_integrate_lg)
  !PRIVATE :: class_adj_tend_lg  ! (called from class_integrate_lg)
  !CLASS--
  !
  ! ISO++
  !PRIVATE :: iso_integrate_gp   ! (called from ch4_physc)
  !PRIVATE :: iso_adj_tend_gp    ! (called from iso_integrate_gp)
  !
!!#D attila +
#ifdef ECHAM5
  !PRIVATE :: iso_integrate_lg   ! (called from ch4_global_end)
  !PRIVATE :: iso_adj_tend_lg    ! (called from iso_integrate_lg)
#endif
!!#D attila -
  ! ISO--

CONTAINS

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE ch4_initialize

    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,     ONLY: find_next_free_unit
    USE messy_main_timer_bi,  ONLY: timer_event_init
    USE messy_ch4,            ONLY: KIE_CH4_13C_OH, &
                                    KIE_CH4_13C_O1D, &
                                    KIE_CH4_13C_CL, &
                                    KIE_CH4_13C_jval, &
                                    KIE_CH4_D1_OH, &
                                    KIE_CH4_D1_O1D, &
                                    KIE_CH4_D1_CL, &
                                    KIE_CH4_D1_jval
    ! SUBROUTINES
    USE messy_ch4,            ONLY: ch4_read_nml_ctrl

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ch4_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

    ! READ CTRL namelist 
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find free I/O unit
       CALL ch4_read_nml_ctrl(status, iou)   ! read CTRL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CTRL namelist',substr)
    END IF
    ! BROADCAST CTRL namelist entries from I/O-PE to ALL OTHER PEs
    ! ISO++
    CALL p_bcast(KIE_CH4_13C_OH, p_io)
    CALL p_bcast(KIE_CH4_13C_O1D, p_io)
    CALL p_bcast(KIE_CH4_13C_CL, p_io)
    CALL p_bcast(KIE_CH4_13C_jval, p_io)
    CALL p_bcast(KIE_CH4_D1_OH, p_io)
    CALL p_bcast(KIE_CH4_D1_O1D, p_io)
    CALL p_bcast(KIE_CH4_D1_CL, p_io)
    CALL p_bcast(KIE_CH4_D1_jval, p_io)

    ! ISO--

    ! READ CPL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find next free I/O unit
       CALL ch4_read_nml_cpl(status, iou)    ! read CPL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CPL namelist',substr)
    END IF
    ! BROADCAST CPL namleist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(i_H2O_feedback, p_io)
    !
    CALL p_bcast(L_GP, p_io)
    CALL p_bcast(c_gp_OH%cha, p_io)
    CALL p_bcast(c_gp_OH%obj, p_io)
    CALL p_bcast(c_gp_O1D%cha, p_io)
    CALL p_bcast(c_gp_O1D%obj, p_io)
    CALL p_bcast(c_gp_Cl%cha, p_io)
    CALL p_bcast(c_gp_Cl%obj, p_io)
    CALL p_bcast(c_gp_jCH4%cha, p_io)
    CALL p_bcast(c_gp_jCH4%obj, p_io)
    ! CLASS++
    CALL p_bcast(i_gp_nclass_emis_age(:), p_io)
    CALL p_bcast(l_gp_adj_tend, p_io)
    CALL p_bcast(i_gp_ageing, p_io)
    CALL p_bcast(r_gp_age_cll, p_io)
    ! CLASS--
    !
    CALL p_bcast(L_LG, p_io)
    CALL p_bcast(c_lg_OH%cha, p_io)
    CALL p_bcast(c_lg_OH%obj, p_io)
    CALL p_bcast(c_lg_O1D%cha, p_io)
    CALL p_bcast(c_lg_O1D%obj, p_io)
    CALL p_bcast(c_lg_Cl%cha, p_io)
    CALL p_bcast(c_lg_Cl%obj, p_io)
    CALL p_bcast(c_lg_jCH4%cha, p_io)
    CALL p_bcast(c_lg_jCH4%obj, p_io)
    ! CLASS++
    CALL p_bcast(i_lg_nclass_emis_age(:), p_io)
    CALL p_bcast(l_lg_adj_tend, p_io)
    CALL p_bcast(i_lg_ageing, p_io)
    CALL p_bcast(r_lg_age_cll, p_io)
    ! CLASS--
    ! ISO++
    CALL p_bcast(l_gp_iso_C, p_io)
    CALL p_bcast(l_gp_iso_H, p_io)
    CALL p_bcast(l_lg_iso_C, p_io)
    CALL p_bcast(l_lg_iso_H, p_io)
    ! ISO--    
    CALL p_bcast(l_ef_re, p_io)  ! op_pj_20160519

    IF ((i_H2O_feedback == 1) .AND. (.NOT. L_GP)) THEN
       CALL error_bi('GP feedback requested, but GP switched off',substr)
    END IF

    IF ((i_H2O_feedback == 2) .AND. (.NOT. L_LG)) THEN
       CALL error_bi('LG feedback requested, but LG switched off',substr)
    END IF

    ! CLASS++
    l_gp_class = (i_gp_nclass_emis_age(1) > 0) .AND. &
                 (i_gp_nclass_emis_age(2) > 0) .AND. L_GP
    l_lg_class = (i_lg_nclass_emis_age(1) > 0) .AND. &
                 (i_lg_nclass_emis_age(2) > 0) .AND. L_LG

    IF (l_gp_class) THEN
       IF (i_gp_nclass_emis_age(2) > NAGECLASSMAX) THEN
          CALL error_bi('Max. number of age classes (GP) exceeded. '//&
               &'Recompile with increased NAGECLASSMAX.',substr)
       END IF
    ENDIF

    IF (l_lg_class) THEN
       IF (i_lg_nclass_emis_age(2) > NAGECLASSMAX) THEN
          CALL error_bi('Max. number of age classes (LG) exceeded. '//&
               &'Recompile with increased NAGECLASSMAX.',substr)
       END IF
    ENDIF

    IF (l_gp_class .OR. l_lg_class) THEN
       CALL timer_event_init(XTIMER_MONTHLY, TIMER_MONTHLY &
            , 'ch4_monthly', 'present')
    END IF
    ! CLASS--

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif

    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

  END SUBROUTINE ch4_initialize
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE ch4_new_tracer

    ! ------------------------------------------------------------------
    ! This subroutine is used to define new tracers. See
    ! http://www.atmos-chem-phys.net/8/1677   (including supplement !)
    ! for full documentation.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
#ifdef ECHAM5 
    USE messy_main_tracer_mem_bi,   ONLY: LGTRSTR
#endif
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    ! MESSy
    USE messy_main_tracer,        ONLY: new_tracer, set_tracer &
                                      , R_molarmass
    USE messy_main_constants_mem, ONLY: MC, MH
    USE messy_main_tools,         ONLY: int2str ! for CLASS

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ch4_new_tracer'
    INTEGER                     :: status
    ! CLASS++
    INTEGER                     :: je, ja
    CHARACTER(LEN=2)            :: se, sa
    INTEGER                     :: idx
    ! CLASS--
    ! ISO++
    INTEGER                         :: ii 
    REAL(DP), DIMENSION(:), POINTER :: iso_molarmass_gp => NULL()
    REAL(DP), DIMENSION(:), POINTER :: iso_molarmass_lg => NULL()
    ! ISO--
    
    CALL start_message_bi(modstr,'TRACER DEFINITION',substr)  ! log-output

    IF (L_GP) THEN
       ! ### define new tracers here
       CALL new_tracer(status, GPTRSTR, 'CH4', modstr, idx=idt_gp_CH4 &
            , subname = 'fx', unit='mol/mol')
       CALL tracer_halt(substr, status)   ! terminate if error
       CALL set_tracer(status, GPTRSTR, idt_gp_CH4 &
            , R_molarmass, r=(MC + 4.0_dp*MH))
       CALL tracer_halt(substr, status)   ! terminate if error

#ifdef MESSYTENDENCY
       CALL mtend_register(my_handle, idt_gp_CH4)
       IF (i_H2O_feedback > 0) &
            CALL mtend_register(my_handle, mtend_id_q)
#endif

       ! CLASS++
       IF (l_gp_class) THEN
          ALLOCATE(idt_gp_class( &
               i_gp_nclass_emis_age(1), i_gp_nclass_emis_age(2) ) )
          idt_gp_class(:,:) = 0
          ALLOCATE(f_gp(i_gp_nclass_emis_age(2)))
          f_gp(:) = 0.0_dp
          emis_gp: DO je=1, i_gp_nclass_emis_age(1)
             CALL int2str(se, je, '0', 'x')
             age_gp: DO ja=1,  i_gp_nclass_emis_age(2)
                CALL int2str(sa, ja, '0', 'x')

                CALL new_tracer(status, GPTRSTR, 'CH4', modstr, idx=idx &
                     , subname = 'fx_e'//se//'_a'//sa &
                     , unit='mol/mol')
                CALL tracer_halt(substr, status)   ! terminate if error
                CALL set_tracer(status, GPTRSTR, idx &
                     , R_molarmass, r=(MC + 4.0_dp*MH))
                CALL tracer_halt(substr, status)   ! terminate if error

#ifdef MESSYTENDENCY
                CALL mtend_register(my_handle, idx)
#endif

                idt_gp_class(je,ja) = idx
               
             END DO age_gp
          END DO emis_gp
       ENDIF
       ! CLASS--
       
       ! ISO++
       ! Set names and molarmasses according to chosen isotopologues
       ! Allocate necessary variables
       IF ((l_gp_iso_C).AND.(l_gp_iso_H)) THEN
          ! carbon and hydrogen isotopologues 
          ALLOCATE(iso_id_gp(4))
          ALLOCATE(iso_molarmass_gp(4))
          iso_id_gp(:) = iso_id(:)
          iso_molarmass_gp(:) = iso_molarmass(:)
       END IF
       IF ((l_gp_iso_C).AND..NOT.(l_gp_iso_H)) THEN
          ! only carbon isotopologues 1, 2
          ALLOCATE(iso_id_gp(2))
          ALLOCATE(iso_molarmass_gp(2))
          iso_id_gp(:) = iso_id(1:2)
          iso_molarmass_gp(:) = iso_molarmass(1:2)
       END IF
       IF (.NOT.(l_gp_iso_C).AND.(l_gp_iso_H)) THEN
          ! only hydrogen isotopologues 3, 4
          ALLOCATE(iso_id_gp(2))
          ALLOCATE(iso_molarmass_gp(2))
          iso_id_gp(:) = iso_id(3:4)
          iso_molarmass_gp(:) = iso_molarmass(3:4)
       END IF

       IF ( l_gp_iso_C .OR. l_gp_iso_H ) THEN

          ALLOCATE(idt_gp_iso(size(iso_id_gp)))
       
          DO ii=1,size(iso_id_gp)
             ! ii = 
             !      1 => 12C
             !      2 => 13C
             !      3 =>  D0
             !      4 =>  D1

             CALL new_tracer(status, GPTRSTR, 'CH4', modstr, &
                  idx=idx, subname = TRIM(iso_name(iso_id_gp(ii))), &
                  unit='mol/mol')
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, GPTRSTR, idx &
                  , R_molarmass, r=iso_molarmass_gp(ii))
             CALL tracer_halt(substr, status)   ! terminate if error
             
             idt_gp_iso(ii) = idx ! save new tracer id
       
             ! HDO++
             IF (iso_id_gp(ii)==iso_D1) THEN
                CALL new_tracer(status, GPTRSTR, 'HDO', modstr, &
                     idx=idx, unit='mol/mol')
                CALL tracer_halt(substr, status)   ! terminate if error
                CALL set_tracer(status, GPTRSTR, idx &
                     , R_molarmass, r=hdo_molarmass)
                CALL tracer_halt(substr, status)   ! terminate if error
                HDO_ex_gp = .TRUE.    ! flag about existence of HDO tracer
                idt_gp_HDO = idx   ! tracer id
                i_gp_CH4_D1 = ii   ! index of CH3D tracer (not its id)         
             END IF
             ! HDO--

#ifdef MESSYTENDENCY
             CALL mtend_register(my_handle,idt_gp_iso(ii))
             ! HDO++
             IF (HDO_ex_gp) &
                  CALL mtend_register(my_handle, idt_gp_HDO)
             ! HDO--
#endif
          END DO
       END IF
       ! ISO--
       
    END IF

!!#D attila +
#ifdef ECHAM5
    IF (L_LG) THEN  
       ! ### define new tracers here
       CALL new_tracer(status, LGTRSTR, 'CH4', modstr, idx=idt_lg_CH4 &
            , subname = 'fx', unit='mol/mol')
       CALL tracer_halt(substr, status)   ! terminate if error
       CALL set_tracer(status, LGTRSTR, idt_lg_CH4 &
            , R_molarmass, r=(MC + 4.0_dp*MH))
       CALL tracer_halt(substr, status)   ! terminate if error

    !qqq+
    ! Note: TENDENCY NOT YET AVAILABLE FOR LG TRACERS
!#ifdef MESSYTENDENCY
!       CALL mtend_register(my_handle, mtend_id_tracer, idt=idt_lg_CH4)
!       IF (i_H2O_feedback > 0) &
!            CALL mtend_register(my_handle, mtend_id_q)
!       IF (idt_lg_H2O > 0) &
!            CALL mtend_register(my_handle, mtend_id_tracer, idt=idt_lg_H2O)
!#endif
    !qqq-

       ! CLASS++
       IF (l_lg_class) THEN
          ALLOCATE(idt_lg_class( &
               i_lg_nclass_emis_age(1), i_lg_nclass_emis_age(2) ) )
          idt_lg_class(:,:) = 0
          ALLOCATE(f_lg(i_lg_nclass_emis_age(2)))
          f_lg(:) = 0.0_dp
          emis_lg: DO je=1, i_lg_nclass_emis_age(1)
             CALL int2str(se, je, '0', 'x')
             age_lg: DO ja=1,  i_lg_nclass_emis_age(2)
                CALL int2str(sa, ja, '0', 'x')

                CALL new_tracer(status, LGTRSTR, 'CH4', modstr, idx=idx &
                     , subname = 'fx_e'//se//'_a'//sa &
                     , unit='mol/mol')
                CALL tracer_halt(substr, status)   ! terminate if error
                CALL set_tracer(status, LGTRSTR, idx &
                     , R_molarmass, r=(MC + 4.0_dp*MH))
                CALL tracer_halt(substr, status)   ! terminate if error

    !qqq+
    ! Note: TENDENCY NOT YET AVAILABLE FOR LG TRACERS
!#ifdef MESSYTENDENCY
!       CALL mtend_register(my_handle, mtend_id_tracer, idt=idx)
!#endif
    !qqq-

                idt_lg_class(je,ja) = idx

             END DO age_lg
          END DO emis_lg
       END IF
       ! CLASS--

       ! ISO++
       ! Set names and molarmasses according to chosen isotopologues
       ! Allocate necessary variables
       IF ((l_lg_iso_C).AND.(l_lg_iso_H)) THEN
          ALLOCATE(iso_id_lg(4))
          ALLOCATE(iso_molarmass_lg(4))
          iso_id_lg(:) = iso_id(:)
          iso_molarmass_lg(:) = iso_molarmass(:)
       END IF
       IF ((l_lg_iso_C).AND..NOT.(l_lg_iso_H)) THEN
          ALLOCATE(iso_id_lg(2))
          ALLOCATE(iso_molarmass_lg(2))
          iso_id_lg(:) = iso_id(1:2)
          iso_molarmass_lg(:) = iso_molarmass(3:4)
       END IF
       IF (.NOT.(l_lg_iso_C).AND.(l_lg_iso_H)) THEN
          ALLOCATE(iso_id_lg(2))
          ALLOCATE(iso_molarmass_lg(2))
          iso_id_lg(:) = iso_id(3:4)
          iso_molarmass_lg(:) = iso_molarmass(3:4)
       END IF

       IF ((l_lg_iso_C).OR.(l_lg_iso_H)) THEN

          ALLOCATE(idt_lg_iso(size(iso_id_lg)))

          DO ii=1,size(iso_id_lg)
             CALL new_tracer(status, LGTRSTR, 'CH4', modstr, &
                  idx=idx, subname = TRIM(iso_name(iso_id_lg(ii))), &
                  unit='mol/mol')
             CALL tracer_halt(substr, status)   ! terminate if error
             CALL set_tracer(status, LGTRSTR, idx &
                  , R_molarmass, r=iso_molarmass_lg(ii))
             CALL tracer_halt(substr, status)   ! terminate if error

             idt_lg_iso(ii) = idx
       
             ! HDO++
             IF (iso_id_lg(ii)==iso_D1) THEN
                CALL new_tracer(status, LGTRSTR, 'HDO', modstr, &
                     idx=idx, unit='mol/mol')
                CALL tracer_halt(substr, status)   ! terminate if error
                CALL set_tracer(status, LGTRSTR, idx &
                     , R_molarmass, r=hdo_molarmass)
                CALL tracer_halt(substr, status)   ! terminate if error
                HDO_ex_lg = .TRUE. ! flag about existence of HDO tracer
                idt_lg_HDO = idx   ! tracer id
                i_lg_CH4_D1 = ii   ! index of CH3D tracer
             END IF
             ! HDO--

             !qqq+
             ! Note: TENDENCY NOT YET AVAILABLE FOR LG TRACERS
!#ifdef MESSYTENDENCY
!             CALL mtend_register(my_handle, mtend_id_tracer, &
!                  idt=idt_lg_iso(ii))
!             ! HDO++
!             IF (HDO_ex_lg) &
!                  CALL mtend_register(my_handle, mtend_id_tracer, &
!                  idt=idt_lg_HDO)
!             ! HDO--
!#endif
          END DO
       END IF
       ! ISO--
    END IF

#endif
!!#D attila -

    ! ISO++
    IF (ASSOCIATED(iso_molarmass_gp)) THEN
       DEALLOCATE(iso_molarmass_gp)
       NULLIFY(iso_molarmass_gp)
    END IF

    IF (ASSOCIATED(iso_molarmass_lg)) THEN
       DEALLOCATE(iso_molarmass_lg)
       NULLIFY(iso_molarmass_lg)
    END IF
    ! ISO--

    CALL end_message_bi(modstr,'TRACER DEFINITION',substr)  ! log-output

  END SUBROUTINE ch4_new_tracer
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE ch4_init_memory

    ! CLASS++
    ! BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
#ifdef ECHAM5
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, LG_ATTILA
#else
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
#endif
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ch4_init_memory'
    INTEGER                     :: status
  
    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    IF (L_GP) THEN
       IF ( (l_gp_class .OR. l_gp_iso_C .OR. l_gp_iso_H) &
            .AND. l_gp_adj_tend ) THEN

          ! define new channel (only once)
          CALL new_channel(status, modstr//'_gp', lrestreq=.FALSE.)
          CALL channel_halt(substr, status)

          IF (l_gp_class) THEN
             CALL new_channel_object(status, modstr//'_gp', 'sadj'&
                  , p3=sadj_gp &
                  , reprid=GP_3D_MID )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_gp' &
                  , 'sadj', 'long_name', c='scaling for tendency adjustment' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_gp' &
                  , 'sadj', 'units', c='-' )
             CALL channel_halt(substr, status)
          END IF

          IF (l_gp_iso_C) THEN
             CALL new_channel_object(status, modstr//'_gp', 'sadj_iso_C' &
                  , p3=sadj_iso_C_gp &
                  , reprid=GP_3D_MID )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_gp' &
                  , 'sadj_iso_C', 'long_name' &
                  , c='scaling for tendency adjustment of carbon isotopologues' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_gp' &
                  , 'sadj_iso_C', 'units', c='-' )
             CALL channel_halt(substr, status)
          END IF
          
          IF (l_gp_iso_H) THEN
             CALL new_channel_object(status, modstr//'_gp', 'sadj_iso_H' &
                  , p3=sadj_iso_H_gp &
                  , reprid=GP_3D_MID )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_gp' &
                  , 'sadj_iso_H', 'long_name' &
                  , c='scaling for tendency adjustment of hydrogen isotopologues' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_gp' &
                  , 'sadj_iso_H', 'units', c='-' )
             CALL channel_halt(substr, status)
          END IF
          
       END IF
    END IF

!!#D attila +
#ifdef ECHAM5
    IF (L_LG) THEN
       IF ( (l_lg_class .OR. l_lg_iso_C .OR. l_lg_iso_H) &
            .AND. l_lg_adj_tend ) THEN

          ! define new channel
          CALL new_channel(status, modstr//'_lg', lrestreq=.FALSE.)
          CALL channel_halt(substr, status)

          IF (l_lg_class) THEN
             CALL new_channel_object(status, modstr//'_lg', 'sadj'&
                  , p1=sadj_lg &
                  , reprid=LG_ATTILA )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_lg' &
                  , 'sadj', 'long_name', c='scaling for tendency adjustment' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_lg' &
                  , 'sadj', 'units', c='-' )
             CALL channel_halt(substr, status)
          END IF
          
          IF (l_lg_iso_C) THEN       
             CALL new_channel_object(status, modstr//'_lg', 'sadj_iso_C' &
                  , p1=sadj_iso_C_lg &
                  , reprid=LG_ATTILA )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_lg' &
                  , 'sadj_iso_C', 'long_name' &
                  , c='scaling for tendency adjustment of carbon isotopologues' )
             
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_lg' &
                  , 'sadj_iso_C', 'units', c='-' )
             CALL channel_halt(substr, status)
          ENDIF

          IF (l_lg_iso_H) THEN
             CALL new_channel_object(status, modstr//'_lg', 'sadj_iso_H' &
                  , p1=sadj_iso_H_lg &
                  , reprid=LG_ATTILA )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_lg' &
                  , 'sadj_iso_H', 'long_name' &
                  , c='scaling for tendency adjustment of hydrogen isotopologues' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_lg' &
                  , 'sadj_iso_H', 'units', c='-' )
             CALL channel_halt(substr, status)
          END IF

       END IF
    END IF
#endif
!!#D attila -

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)
    ! CLASS--

  END SUBROUTINE ch4_init_memory
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE ch4_init_coupling

    ! ------------------------------------------------------------------
    ! This soubroutine is used to set pointers
    ! (channel objects and/or tracers) for coupling to the 
    ! basemodel and to other submodels.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: get_channel_object &
                                         , set_channel_object_restreq
#ifdef ECHAM5
    USE messy_main_tracer_mem_bi,    ONLY: LGTRSTR
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_tracer,           ONLY: get_tracer, TR_NEXIST
#endif
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ch4_init_coupling'
    INTEGER                     :: status

    CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output

    IF (L_GP) THEN

       CALL get_channel_object(status &
            , TRIM(c_gp_OH%cha), TRIM(c_gp_OH%obj), p3=ptr_gp_OH)
       CALL channel_halt(substr, status)

       CALL warning_bi( &
            'restart flag activated for object '//&
            &TRIM(c_gp_OH%obj)//' in channel '//TRIM(c_gp_OH%cha), substr)
       CALL set_channel_object_restreq(status, &
            TRIM(c_gp_OH%cha), TRIM(c_gp_OH%obj) )
       CALL channel_halt(substr, status)
       
       ! ---

       CALL get_channel_object(status &
            , TRIM(c_gp_O1D%cha), TRIM(c_gp_O1D%obj), p3=ptr_gp_O1D)
       CALL channel_halt(substr, status)

       CALL warning_bi( &
            'restart flag activated for object '//&
            &TRIM(c_gp_O1D%obj)//' in channel '//TRIM(c_gp_O1D%cha), substr)
       CALL set_channel_object_restreq(status, &
            TRIM(c_gp_O1D%cha), TRIM(c_gp_O1D%obj) )
       CALL channel_halt(substr, status)

       ! ---

       CALL get_channel_object(status &
            , TRIM(c_gp_Cl%cha), TRIM(c_gp_Cl%obj), p3=ptr_gp_Cl)
       CALL channel_halt(substr, status)

       CALL warning_bi( &
            'restart flag activated for object '//&
            &TRIM(c_gp_Cl%obj)//' in channel '//TRIM(c_gp_Cl%cha), substr)
       CALL set_channel_object_restreq(status, &
            TRIM(c_gp_Cl%cha), TRIM(c_gp_Cl%obj) )
       CALL channel_halt(substr, status)
       
       ! ---

       CALL get_channel_object(status &
            , TRIM(c_gp_jCH4%cha), TRIM(c_gp_jCH4%obj), p3=ptr_gp_jCH4)
       CALL channel_halt(substr, status)

       CALL warning_bi( &
            'restart flag activated for object '//&
            &TRIM(c_gp_jCH4%obj)//' in channel '//TRIM(c_gp_jCH4%cha), substr)
       CALL set_channel_object_restreq(status, &
            TRIM(c_gp_jCH4%cha), TRIM(c_gp_jCH4%obj) )
       CALL channel_halt(substr, status)
    END IF
       
!!#D attila +
#ifdef ECHAM5
    IF (L_LG) THEN
       CALL get_channel_object(status &
            , 'attila', 'PPRESS', p1=press_3d_lg)
       CALL channel_halt(substr, status)

       CALL get_tracer(status, LGTRSTR, 'H2O', idx=idt_lg_H2O)
       IF (status == TR_NEXIST) THEN
          CALL warning_bi('LAGRANGIAN H2O TRACER NOT PRESENT!', substr)
          idt_lg_H2O = 0
       ELSE
          CALL tracer_halt(substr, status) 
       END IF

       ! ---

       CALL get_channel_object(status &
            , TRIM(c_lg_OH%cha), TRIM(c_lg_OH%obj), p1=ptr_lg_OH)
       CALL channel_halt(substr, status)

       CALL warning_bi( &
            'restart flag activated for object '//&
            &TRIM(c_lg_OH%obj)//' in channel '//TRIM(c_lg_OH%cha), substr)
       CALL set_channel_object_restreq(status, &
            TRIM(c_lg_OH%cha), TRIM(c_lg_OH%obj) )
       CALL channel_halt(substr, status)

       ! ---
       
       CALL get_channel_object(status &
            , TRIM(c_lg_O1D%cha), TRIM(c_lg_O1D%obj), p1=ptr_lg_O1D)
       CALL channel_halt(substr, status)

       CALL warning_bi( &
            'restart flag activated for object '//&
            &TRIM(c_lg_O1D%obj)//' in channel '//TRIM(c_lg_O1D%cha), substr)
       CALL set_channel_object_restreq(status, &
            TRIM(c_lg_O1D%cha), TRIM(c_lg_O1D%obj) )
       CALL channel_halt(substr, status)
       
       ! ---

       CALL get_channel_object(status &
            , TRIM(c_lg_Cl%cha), TRIM(c_lg_Cl%obj), p1=ptr_lg_Cl)
       CALL channel_halt(substr, status)

       CALL warning_bi( &
            'restart flag activated for object '//&
            &TRIM(c_lg_Cl%obj)//' in channel '//TRIM(c_lg_Cl%cha), substr)
       CALL set_channel_object_restreq(status, &
            TRIM(c_lg_Cl%cha), TRIM(c_lg_Cl%obj) )
       CALL channel_halt(substr, status)
       
       ! ---

       CALL get_channel_object(status &
            , TRIM(c_lg_jCH4%cha), TRIM(c_lg_jCH4%obj), p1=ptr_lg_jCH4)
       CALL channel_halt(substr, status)

       CALL warning_bi( &
            'restart flag activated for object '//&
            &TRIM(c_lg_jCH4%obj)//' in channel '//TRIM(c_lg_jCH4%cha), substr)
       CALL set_channel_object_restreq(status, &
            TRIM(c_lg_jCH4%cha), TRIM(c_lg_jCH4%obj) )
       CALL channel_halt(substr, status)

       ! ---

    END IF
#endif
!!#D attila -

    CALL end_message_bi(modstr,'COUPLING',substr)  ! log-output

  END SUBROUTINE ch4_init_coupling
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE ch4_global_start

    ! CLASS++
    USE messy_main_timer_bi,      ONLY: event_state
    USE messy_main_timer,         ONLY: current_date, time_step_len

    IMPLICIT NONE
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ch4_global_start'

    IF (l_gp_class) THEN
       SELECT CASE(i_gp_ageing)
       CASE(0,2)
          ! case 0 is NOT recommended (see note in class_age_move_gp)
          ! case 2 cannot be used with MESSYTENDENCY
          l_gp_mv_agecl_now = event_state(XTIMER_MONTHLY, current_date)
          f_gp(:) = 1.0_dp
       CASE(1)
          l_gp_mv_agecl_now = .TRUE.
          ! convert from days to seconds
          f_gp = time_step_len / &
               (r_gp_age_cll(1:i_gp_nclass_emis_age(2))*86400.0_dp)
       CASE DEFAULT
          ! other cases might require additional changes in class_age_move_gp
          CALL error_bi('unknown ageing method (GP)',substr)
       END SELECT
    END IF

    IF (l_lg_class) THEN
       SELECT CASE(i_lg_ageing)
       CASE(0,2)
          ! case 0 is NOT recommended (see note in class_age_move_gp)
          ! case 2 cannot be used with MESSYTENDENCY
          l_lg_mv_agecl_now = event_state(XTIMER_MONTHLY, current_date)
          f_lg = 1.0_dp
       CASE(1)
          ! convert from days to seconds
          f_lg = time_step_len / &
               (r_lg_age_cll(1:i_lg_nclass_emis_age(2))*86400.0_dp)    
       CASE DEFAULT
          ! other cases might require additional changes in class_age_move_lg
          CALL error_bi('unknown ageing method (LG)',substr)
       END SELECT
    END IF

    ! CLASS--

  END SUBROUTINE ch4_global_start
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE ch4_vdiff

  END SUBROUTINE ch4_vdiff
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE ch4_physc

    USE messy_main_constants_mem,   ONLY: M_air, M_H2O
    USE messy_main_timer,           ONLY: time_step_len
    USE messy_main_grid_def_mem_bi, ONLY: kproma, nlev, jrow
    USE messy_main_data_bi,         ONLY: press_3d, tm1_3d, tte_3d  &
                                        , qm1_3d, qte_3d
#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, pxtm1 => qxtm1
#endif

! Necessary: iso_ratio_* > ratio of isotope tracer in overall budget
!            eps         > relative tolerance for subtracers (e.g. 0.001*CH4_te)

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*),       PARAMETER :: substr = 'ch4_physc'
    REAL(DP),               PARAMETER :: scvmr = M_air/M_H2O
    INTEGER                           :: idt
    REAL(DP), DIMENSION(:,:), POINTER :: temp     => NULL()
    REAL(DP), DIMENSION(:,:), POINTER :: spechum  => NULL()
    REAL(DP), DIMENSION(:,:), POINTER :: CH4      => NULL()
    REAL(DP), DIMENSION(:,:), POINTER :: CH4_te   => NULL() 
    REAL(DP), DIMENSION(:,:), POINTER :: q_te     => NULL() 

    REAL(DP), DIMENSION(:,:), POINTER :: h2o_bc   => NULL() 
    REAL(DP), DIMENSION(:,:), POINTER :: h2o_ac   => NULL() 
    REAL(DP), DIMENSION(:,:), POINTER :: q_ac     => NULL() 

    IF (.NOT. L_GP) RETURN

    ! LOCAL SPACE
    ALLOCATE(temp(kproma, nlev))
    ALLOCATE(spechum(kproma, nlev))
    ALLOCATE(CH4(kproma, nlev))
    ALLOCATE(CH4_te(kproma, nlev))
    IF (i_H2O_feedback == 1) THEN
       ALLOCATE(q_te(kproma, nlev))
       ALLOCATE(h2o_bc(kproma, nlev))
       ALLOCATE(h2o_ac(kproma, nlev))
       ALLOCATE(q_ac(kproma, nlev))
    END IF

    ! CALUCLATE CORRECT START VALUES
#ifndef MESSYTENDENCY
    temp(:,:)  = tm1_3d(_RI_XYZ__(1:kproma,jrow,:)) &
         + tte_3d(_RI_XYZ__(1:kproma,jrow,:)) * time_step_len

    spechum(:,:)  = qm1_3d(_RI_XYZ__(1:kproma,jrow,:)) &
         + qte_3d(_RI_XYZ__(1:kproma,jrow,:)) * time_step_len

    idt = idt_gp_CH4
    CH4(:,:) = pxtm1(_RI_X_ZN_(1:kproma,:,idt)) &
         + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len

#else
    CALL mtend_get_start_l(mtend_id_t, v0=temp(:,:))
    CALL mtend_get_start_l(mtend_id_q, v0=spechum(:,:))
    CALL mtend_get_start_l(idt_gp_CH4, v0=CH4(:,:))
#endif
    CH4_te(:,:) = 0.0_dp 

    CALL ch4_integrate(&
         CH4_te(1:kproma,:)            &
         , CH4(1:kproma,:)             &
         , ptr_gp_OH(_RI_XYZ__(1:kproma,jrow,:))    &
         , ptr_gp_O1D(_RI_XYZ__(1:kproma,jrow,:))   &
         , ptr_gp_Cl(_RI_XYZ__(1:kproma,jrow,:))    &
         , ptr_gp_jCH4(_RI_XYZ__(1:kproma,jrow,:))  &
         , temp(1:kproma,:)            &
         , press_3d(_RI_XYZ__(1:kproma,jrow,:))  &
         , spechum(1:kproma,:)         &
         )

    IF (i_H2O_feedback == 1) THEN
       !
       !
       !        dq    d      H2O           scvmr * d(H2O)/dt
       ! qte = ---- = -- (-----------) = ---------------------
       !        dt    dt  scvmr + H2O      (scvmr + H2O)^2
       !
       ! d(H2O)/dt = -2*d(CH4)/dt ! ONE METHANE MOLECULE REACTS INTO 
       !                          ! 2 WATER MOLECULES
       !                 q
       ! H2O = scvmr * -----
       !               1 - q
       !
       h2o_bc(1:kproma,:) = scvmr * spechum(1:kproma,:) / &
            (1.0_dp - spechum(1:kproma,:))
       h2o_ac(1:kproma,:) = h2o_bc(1:kproma,:) - 2._dp*CH4_te(1:kproma,:)
       q_ac(1:kproma,:)   = h2o_ac(1:kproma,:) / (scvmr + h2o_ac(1:kproma,:))
       q_te(1:kproma,:)   = ( q_ac(1:kproma,:) - spechum(1:kproma,:) ) / &
            time_step_len
    END IF

    ! Update tracers
    ! Check for consistency of Subtracers
#ifndef MESSYTENDENCY
    idt = idt_gp_CH4 ! main tracer
    pxtte(_RI_X_ZN_(1:kproma,:,idt)) = pxtte(_RI_X_ZN_(1:kproma,:,idt)) + CH4_te(1:kproma,:)

    ! H2O tracer
    IF (i_H2O_feedback == 1) THEN
       qte_3d(_RI_XYZ__(1:kproma,jrow,:)) = qte_3d(_RI_XYZ__(1:kproma,jrow,:)) + q_te(1:kproma,:)
    END IF
#else
    ! main tracer
    CALL mtend_add_l(my_handle, idt_gp_CH4, px=CH4_te)

    ! H2O tracer
    IF (i_H2O_feedback == 1) THEN
       CALL mtend_add_l(my_handle, mtend_id_q, px=q_te)
    END IF
#endif

    ! CLASS++
    IF (l_gp_class) &
         CALL class_integrate_gp( &
         temp(1:kproma,:), press_3d(_RI_XYZ__(1:kproma,jrow,:)), spechum(1:kproma,:) )
    ! CLASS--

    ! ISO++
    IF (l_gp_iso_C.OR.l_gp_iso_H) &
         CALL iso_integrate_gp( &
         temp(1:kproma,:), press_3d(_RI_XYZ__(1:kproma,jrow,:)), spechum(1:kproma,:), CH4_te(1:kproma, :) )
    ! ISO--

    ! CLEAN MEMORY
    DEALLOCATE(temp)    ; NULLIFY(temp)
    DEALLOCATE(spechum) ; NULLIFY(spechum)
    DEALLOCATE(CH4)     ; NULLIFY(CH4)
    DEALLOCATE(CH4_te)  ; NULLIFY(CH4_te)

    IF (i_H2O_feedback == 1) THEN
       DEALLOCATE(q_te) ; NULLIFY(q_te)
       DEALLOCATE(h2o_bc) ; NULLIFY(h2o_bc)
       DEALLOCATE(h2o_ac) ; NULLIFY(h2o_ac)
       DEALLOCATE(q_ac) ; NULLIFY(q_ac)
    END IF
    
  END SUBROUTINE ch4_physc
  ! =========================================================================

  ! === PRIVATE =============================================================
  SUBROUTINE class_integrate_gp(temp, press, spechum)

    USE messy_main_timer,           ONLY: time_step_len
    USE messy_main_grid_def_mem_bi, ONLY: kproma, nlev, jrow
#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi,   ONLY: pxtte => qxtte, pxtm1 => qxtm1
#endif
    USE messy_main_tools,         ONLY: PTR_2D_ARRAY

    IMPLICIT NONE

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: temp
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: press
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: spechum

    ! LOCAL
    INTEGER                                     :: je, ja 
    INTEGER                                     :: idt
    TYPE(PTR_2D_ARRAY), DIMENSION(:,:), POINTER :: CH4c    => NULL()
    TYPE(PTR_2D_ARRAY), DIMENSION(:,:), POINTER :: CH4c_te => NULL()

    ! allocate space for tracers and tendencies
    ! (need to be stored individually for adjustment)
    ALLOCATE(CH4c(i_gp_nclass_emis_age(1), i_gp_nclass_emis_age(2)))
    ALLOCATE(CH4c_te(i_gp_nclass_emis_age(1), i_gp_nclass_emis_age(2)))

    emis_gp: DO je=1, i_gp_nclass_emis_age(1)
       age_gp: DO ja=1,  i_gp_nclass_emis_age(2)

          ! allocate space for tracers and tendencies
          ! (need to be stored individually for adjustment)
          ALLOCATE(CH4c(je,ja)%ptr(kproma,nlev))
          CH4c(je,ja)%ptr(:,:) = 0.0_dp
          ALLOCATE(CH4c_te(je,ja)%ptr(kproma,nlev))
          CH4c_te(je,ja)%ptr(:,:) = 0.0_dp

          ! calculate individual start values
#ifndef MESSYTENDENCY
          idt = idt_gp_class(je,ja)
          CH4c(je,ja)%ptr(:,:) = pxtm1(_RI_X_ZN_(1:kproma,:,idt)) &
               + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
#else
          CALL mtend_get_start_l(idt_gp_class(je,ja), v0=CH4c(je,ja)%ptr(:,:))
#endif

          ! integrate chemistry and store individual tendencies
          CALL ch4_integrate(&
               CH4c_te(je,ja)%ptr(1:kproma,:)   &
               , CH4c(je,ja)%ptr(1:kproma,:)    &
               , ptr_gp_OH(_RI_XYZ__(1:kproma,jrow,:))    &
               , ptr_gp_O1D(_RI_XYZ__(1:kproma,jrow,:))   &
               , ptr_gp_Cl(_RI_XYZ__(1:kproma,jrow,:))    &
               , ptr_gp_jCH4(_RI_XYZ__(1:kproma,jrow,:))  &
               , temp(:,:)                      &
               , press(:,:)                     &
               , spechum(:,:)                   &
               )

       END DO age_gp
    END DO emis_gp

    CALL class_age_move_gp(CH4c, CH4c_te)

    IF (l_gp_adj_tend) CALL class_adj_tend_gp(CH4c, CH4c_te)

    ! ADD TRACER TENDENCIES
    DO je=1, i_gp_nclass_emis_age(1)
       DO ja=1,  i_gp_nclass_emis_age(2)
#ifndef MESSYTENDENCY
          idt = idt_gp_class(je,ja)
          pxtte(_RI_X_ZN_(1:kproma,:,idt)) = &
               pxtte(_RI_X_ZN_(1:kproma,:,idt)) + CH4c_te(je,ja)%ptr(:,:)
#else
          CALL mtend_add_l(my_handle, idt_gp_class(je,ja) &
               , px=CH4c_te(je,ja)%ptr)
#endif

          ! not longer needed, clean memory
          DEALLOCATE(CH4c(je,ja)%ptr)   ; NULLIFY(CH4c(je,ja)%ptr)
          DEALLOCATE(CH4c_te(je,ja)%ptr); NULLIFY(CH4c_te(je,ja)%ptr)
       END DO
    END DO

    ! CLEAN MEMORY
    DEALLOCATE(CH4c)      ; NULLIFY(CH4c)
    DEALLOCATE(CH4c_te)   ; NULLIFY(CH4c_te)

  END SUBROUTINE class_integrate_gp
  ! =========================================================================

  ! === PRIVATE =============================================================
  SUBROUTINE class_age_move_gp(CH4c, CH4c_te)

    ! UPDATE AGE CLASSES

    USE messy_main_tools,         ONLY: PTR_2D_ARRAY
    USE messy_main_timer,         ONLY: time_step_len

    USE messy_main_tracer_mem_bi, ONLY: qxt, qxtm1, qxtf, qxtte
    USE messy_main_data_bi,         ONLY: eps
    USE messy_main_grid_def_mem_bi, ONLY: kproma

    IMPLICIT NONE

    ! I/O
    TYPE(PTR_2D_ARRAY), DIMENSION(:,:), POINTER :: CH4c
    TYPE(PTR_2D_ARRAY), DIMENSION(:,:), POINTER :: CH4c_te

    ! LOCAL
    INTEGER  :: ja, je ! age and emission class counter
    REAL(DP) :: sro    ! scale for 'remove own'
    INTEGER  :: idt, idt2   ! op_pj_20140717

    !
    ! EMISSIONS GO TO AGE CLASS eXX_a01 ... (youngest!)
    !
    ! A' = A + M*A
    !
    !                  = TENDENCY * dt
    !                |---------------------------------------------|
    ! A'(1) = A(1) +  -1*A(1)                                      |
    ! A'(2) = A(2) +   1*A(1) -1*A(2)                              |
    ! A'(3) = A(3) +           1*A(2) -1*A(3)                      | *dt
    ! A'(4) = A(4) +                   1*A(3) -1*A(4)              |
    ! ...                                                          |
    ! A'(n) = A(n) +                           + 1*A(n-1) -0*A(n)  |
    !
    ! INSTANTANEOUS "STEP": TENDENCY = M*A / dt
    !
    ! => M =  [  -1  0  0  0  0 ... ]    A(1)
    !         [   1 -1  0  0  0 ... ]    A(2)
    !         [   0  1 -1  0  0 ... ]    A(3)
    !         [   0  0  1 -1  0 ... ] *  A(4)
    !         [   ...           ... ]    ...
    !         [   .. 0  0  0  1   0 ]    A(n)
    !
    ! Note: according to the Leapfrog time stepping with
    !       Asselin-filter, this might cause numerical
    !       oszillations with negative values etc. ...
    !
    !
    ! QUASI CONTINUOUS: f = dt/DELTA-t (DELTA-t = age class time span)
    !
    ! A(1) = A(1) + -f*A(1)
    ! A(2) = A(2) +  f*A(1) -f*A(2)
    ! A(3) = A(3) +          f*A(2) -f*A(3)
    ! A(4) = A(4) +                  f*A(3) -f*A(4)
    ! ...
    ! A(n) = A(n) +                          f*A(n-1) -0*A(n)
    !
    ! => M' = [ -f  0  0  0  0 ... ] = f*M
    !         [  f -f  0  0  0 ... ]
    !         [  0  f -f  0  0 ... ]
    !         [  0  0  f -f  0 ... ]
    !         [  ...           ... ]
    !         [  ... 0 0  0  f   0 ]
    !

    IF (.NOT. (l_gp_mv_agecl_now .AND. (i_gp_nclass_emis_age(2) > 1))) RETURN
    
    ! ELSE end of month is reached and age classes are more than one
    SELECT CASE(i_gp_ageing)

    CASE(0,1)

       DO ja = i_gp_nclass_emis_age(2), 2, -1
          IF (ja == i_gp_nclass_emis_age(2)) THEN
             sro = 0.0_dp
          ELSE
             sro = -1.0_dp
          END IF
          DO je=1, i_gp_nclass_emis_age(1)
             ! tendency to remove entire own class  (not for oldest!)
             ! + tendency to add entire younger class
             CH4c_te(je,ja)%ptr(:,:) = &
                  f_gp(ja) * ( &
                  sro * ( ( CH4c(je,ja)%ptr(:,:) &
                  + CH4c_te(je,ja)%ptr(:,:) * time_step_len )/time_step_len ) )&
                  + f_gp(ja-1) * ( ( CH4c(je,ja-1)%ptr(:,:) &
                  + CH4c_te(je,ja-1)%ptr(:,:) * time_step_len )/time_step_len )
          END DO
       END DO
       ja = 1
       DO je=1, i_gp_nclass_emis_age(1)
          CH4c_te(je,ja)%ptr(:,:) = &
               f_gp(ja) * ( &
               -1.0_dp * ( ( CH4c(je,ja)%ptr(:,:) &
               + CH4c_te(je,ja)%ptr(:,:) * time_step_len )/time_step_len ) &
               )
       END DO

    CASE(2)
       ! RESET TRACERS COMPLETELY
       DO ja = i_gp_nclass_emis_age(2), 2, -1      ! loop over age class
          IF (ja == i_gp_nclass_emis_age(2)) THEN  ! last age class, no removal
             sro = 0.0_dp
          ELSE
             sro = -1.0_dp ! removal of former age class
          END IF
          DO je=1, i_gp_nclass_emis_age(1)  ! loop over emissions
             idt  = idt_gp_class(je,ja)     ! tracer identifiers
             idt2 = idt_gp_class(je,ja-1)   ! identifier of former ageclass

             ! 1) everything what's been in the class before
             ! 2) removal of everything what's been in the class before,
             !    if it is not the last class
             ! 3) add everything what's been in the younger class
             IF (ASSOCIATED(qxt)) & 
                  qxt(_RI_X_ZN_(1:kproma,:,idt)) = qxt(_RI_X_ZN_(1:kproma,:,idt)) + & ! 1)
                  sro * qxt(_RI_X_ZN_(1:kproma,:,idt)) + &                   ! 2)
                  qxt(_RI_X_ZN_(1:kproma,:,idt2))                            ! 3)
             IF (ASSOCIATED(qxtm1)) &
                  qxtm1(_RI_X_ZN_(1:kproma,:,idt)) = (1.0_dp - eps) *  &
                  qxt(_RI_X_ZN_(1:kproma,:,idt))
             IF (ASSOCIATED(qxtf)) &
                  qxtf(_RI_X_ZN_(1:kproma,:,idt)) = qxt(_RI_X_ZN_(1:kproma,:,idt)) + & 
                  eps*(qxtm1(_RI_X_ZN_(1:kproma,:,idt)) - 2.0_dp*qxt(_RI_X_ZN_(1:kproma,:,idt)))
             IF (ASSOCIATED(qxtte)) &        ! tendency
                  qxtte(_RI_X_ZN_(1:kproma,:,idt)) = qxtte(_RI_X_ZN_(1:kproma,:,idt)) + & 
                  sro * qxtte(_RI_X_ZN_(1:kproma,:,idt)) + &
                  qxtte(_RI_X_ZN_(1:kproma,:,idt2))
          END DO
       END DO
       ja = 1
       DO je=1, i_gp_nclass_emis_age(1)
          ! set mixing ratio in first age class for all emission 
          ! classes to zero
          idt  = idt_gp_class(je,ja)
          IF (ASSOCIATED(qxt))   qxt(_RI_X_ZN_(1:kproma,:,idt))   = 0.0_dp ! mix. ratio
          IF (ASSOCIATED(qxtm1)) qxtm1(_RI_X_ZN_(1:kproma,:,idt)) = 0.0_dp ! at t-1
          IF (ASSOCIATED(qxtf))  qxtf(_RI_X_ZN_(1:kproma,:,idt))  = 0.0_dp ! filtered
          IF (ASSOCIATED(qxtte)) qxtte(_RI_X_ZN_(1:kproma,:,idt)) = 0.0_dp ! tendency
       END DO

    END SELECT
    
  END SUBROUTINE class_age_move_gp
  ! =========================================================================

  ! === PRIVATE =============================================================
  SUBROUTINE class_adj_tend_gp(CH4c, CH4c_te)

    USE messy_main_grid_def_mem_bi, ONLY: kproma, nlev, jrow
    USE messy_main_timer,           ONLY: time_step_len
    USE messy_main_tracer_mem_bi,   ONLY: pxtte => qxtte, pxtm1 => qxtm1
    USE messy_main_tools,           ONLY: PTR_2D_ARRAY

    IMPLICIT NONE

    ! I/O
    TYPE(PTR_2D_ARRAY), DIMENSION(:,:), POINTER :: CH4c
    TYPE(PTR_2D_ARRAY), DIMENSION(:,:), POINTER :: CH4c_te
    
    ! LOCAL
    INTEGER                           :: je, ja
    INTEGER                           :: idt
    REAL(DP), DIMENSION(:,:), POINTER :: CH4_m1_sum   => NULL() 
    REAL(DP), DIMENSION(:,:), POINTER :: CH4_te_sum   => NULL() 
    REAL(DP), DIMENSION(:,:), POINTER :: tmp => NULL()

    ! --- ADJUST TENDENCIES TO MASTER TRACER .......................
    ! --- 1) MEMORY
    ALLOCATE(CH4_m1_sum(kproma, nlev))
    CH4_m1_sum(:,:) = 0.0_dp
    ALLOCATE(CH4_te_sum(kproma, nlev))
    CH4_te_sum(:,:) = 0.0_dp
    ALLOCATE(tmp(kproma, nlev))
    tmp(:,:) = 0.0_dp
    !
    sadj_gp(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    !
    ! --- 2) calculate SUM
    DO je=1, i_gp_nclass_emis_age(1)
       DO ja=1,  i_gp_nclass_emis_age(2)
          idt = idt_gp_class(je,ja)
          CH4_m1_sum(:,:) = CH4_m1_sum(:,:) + pxtm1(_RI_X_ZN_(1:kproma,:,idt))
          CH4_te_sum(:,:) = CH4_te_sum(:,:) + pxtte(_RI_X_ZN_(1:kproma,:,idt)) + &
               CH4c_te(je,ja)%ptr(:,:)   ! not yet added !!!
       END DO
    END DO
    !
    ! --- 3) calculate scaling (weight) factor
    ! Note: 'CH4' tendency has already been added above to master tracer ...
    idt = idt_gp_CH4 ! master tracer
    CALL sca_tend(pxtm1(_RI_X_ZN_(1:kproma,:,idt)), pxtte(_RI_X_ZN_(1:kproma,:,idt)) &
         , CH4_m1_sum(:,:), CH4_te_sum(:,:), time_step_len &
         , sadj_gp(_RI_XYZ__(1:kproma,jrow,:)))
    !
    ! --- 4) apply scaling / add additional tendency for adjustment
    DO je=1, i_gp_nclass_emis_age(1)
       DO ja=1,  i_gp_nclass_emis_age(2)
          idt = idt_gp_class(je,ja)
          CALL adj_tend(pxtm1(_RI_X_ZN_(1:kproma,:,idt))  &
               , pxtte(_RI_X_ZN_(1:kproma,:,idt))         &
               , sadj_gp(_RI_XYZ__(1:kproma,jrow,:))       &
               , time_step_len                   &
               , tmp(:,:)                        &
               )
          CH4c_te(je,ja)%ptr(:,:) = CH4c_te(je,ja)%ptr(:,:) + tmp(:,:)
       END DO
    END DO
    ! --- END ADJUST TENDENCIES TO MASTER TRACER ....................

    ! CLEAN MEMORY
    DEALLOCATE(tmp)       ; NULLIFY(tmp)
    DEALLOCATE(CH4_m1_sum); NULLIFY(CH4_m1_sum)
    DEALLOCATE(CH4_te_sum); NULLIFY(CH4_te_sum)

  END SUBROUTINE class_adj_tend_gp
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE ch4_global_end

!!#D attila +
#ifdef ECHAM5
    USE messy_main_constants_mem,   ONLY: M_air, M_H2O
    USE messy_main_timer,           ONLY: time_step_len
    USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev, ngpblks
    USE messy_main_data_bi,         ONLY: tm1_3d, tte_3d  &
                                        , qm1_3d, qte_3d
    USE messy_main_tracer_mem_bi,   ONLY: qxtte_a, qxtm1_a &
                                         , NCELL
    USE messy_attila_tools_e5,      ONLY: gp2lg_e5, lg2gp_e5, LG2GP_SUM

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*),         PARAMETER :: substr = 'ch4_global_end'
    REAL(DP),                 PARAMETER :: scvmr = M_air/M_H2O
    INTEGER                             :: idt
    REAL(DP), DIMENSION(:,:,:), POINTER :: temp        => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: spechum     => NULL()
    REAL(DP), DIMENSION(:),     POINTER :: temp_lg     => NULL()
    REAL(DP), DIMENSION(:),     POINTER :: spechum_lg  => NULL()
    !
    REAL(DP), DIMENSION(:),     POINTER :: CH4_lg    => NULL()
    REAL(DP), DIMENSION(:),     POINTER :: CH4_te_lg => NULL() 
    REAL(DP), DIMENSION(:),     POINTER :: q_te_lg   => NULL() 
    REAL(DP), DIMENSION(:,:,:), POINTER :: q_te      => NULL() 
    INTEGER :: jrow

    IF (.NOT. L_LG) RETURN

    ! LOCAL SPACE
    ALLOCATE(temp(nproma, nlev, ngpblks))
    ALLOCATE(spechum(nproma, nlev, ngpblks))
    !
    ALLOCATE(temp_lg(NCELL))
    ALLOCATE(spechum_lg(NCELL))
    ALLOCATE(CH4_lg(NCELL))
    ALLOCATE(CH4_te_lg(NCELL))
    !
    IF (i_H2O_feedback == 2) THEN
       ALLOCATE(q_te(nproma, nlev, ngpblks))
       ALLOCATE(q_te_lg(NCELL))
    END IF

    ! CALUCLATE CORRECT START VALUES
#ifndef MESSYTENDENCY
    DO jrow=1, ngpblks
       temp(:,:,jrow)    = tm1_3d(_RI_XYZ__(:,jrow,:)) + tte_3d(_RI_XYZ__(:,jrow,:)) * time_step_len
       spechum(:,:,jrow) = qm1_3d(_RI_XYZ__(:,jrow,:)) + qte_3d(_RI_XYZ__(:,jrow,:)) * time_step_len
    END DO

    idt = idt_lg_CH4
    CH4_lg(:) = qxtm1_a(:,idt) + qxtte_a(:,idt) * time_step_len       
#else
    CALL mtend_get_start_g(mtend_id_t, v0=temp(:,:,:))
    CALL mtend_get_start_g(mtend_id_q, v0=spechum(:,:,:))
    ! qqq+ NOTE: TENDENCY NOT YET IMPLEMENTED FOR LG TRACERS ...
    !!$CALL mtend_get_start_g(mtend_id_tracer, v0=CH4(:), idt=idt_lg_CH4)
    idt = idt_lg_CH4
    CH4_lg(:) = qxtm1_a(:,idt) + qxtte_a(:,idt) * time_step_len       
    ! qqq-
#endif
    CH4_te_lg(:) = 0.0_dp 

    ! TRANSFORM TO LG
    CALL gp2lg_e5(temp, temp_lg)
    CALL gp2lg_e5(spechum, spechum_lg)

    CALL ch4_integrate(    &
         CH4_te_lg(:)      &
         , CH4_lg(:)       &
         , ptr_lg_OH(:)    &
         , ptr_lg_O1D(:)   &
         , ptr_lg_Cl(:)    &
         , ptr_lg_jCH4(:)  &
         , temp_lg(:)      &
         , press_3d_lg(:)  &
         , spechum_lg(:)   &
         )

    IF (i_H2O_feedback == 2) THEN
       ! ONE METHANE MOLECULE REACTS INTO 2 WATER MOLECULES
       q_te_lg(:) = (-2._dp*CH4_te_lg(:)) &
            / (scvmr * (1.0_dp / (1.0_dp - spechum_lg(:)))**2) 
       CALL lg2gp_e5(q_te_lg, q_te, LG2GP_SUM, fill_value=0.0_dp)
    END IF

#ifndef MESSYTENDENCY
    idt = idt_lg_CH4
    qxtte_a(:,idt) = qxtte_a(:,idt) + CH4_te_lg(:)
    IF (i_H2O_feedback == 2) THEN
       DO jrow=1, ngpblks
          qte_3d(_RI_XYZ__(:,jrow,:)) = qte_3d(_RI_XYZ__(:,jrow,:)) + q_te(:,:,jrow)
       END DO
    END IF
    ! MODIFY LAGRANGIAN H2O TRACER, IF PRESENT
    IF (idt_lg_H2O > 0) THEN
       ! ONE METHANE MOLECULE REACTS INTO 2 WATER MOLECULES
       qxtte_a(:,idt_lg_H2O) = qxtte_a(:,idt_lg_H2O) - 2._dp*CH4_te_lg(:)
    END IF
#else
! qqq+
! NOTE: TENDENCY NOT YET IMPLEMENTED FOR LG TRACER 
    !!$ CALL mtend_add_g(my_handle, mtend_id_tracer&
    !!$                , px=CH4_te_lg, idt=idt_lg_CH4)
    idt = idt_lg_CH4
    qxtte_a(:,idt) = qxtte_a(:,idt) + CH4_te_lg(:)
! qqq-

    ! MODIFY LAGRANGIAN H2O TRACER, IF PRESENT
    IF (idt_lg_H2O > 0) THEN
! qqq+
! NOTE: TENDENCY NOT YET IMPLEMENTED FOR LG TRACER 
    !!$ CALL mtend_add_g(my_handle, mtend_id_tracer &
    !!$                , px=-2._dp*CH4_te_lg, idt=idt_lg_H2O)
       ! ONE METHANE MOLECULE REACTS INTO 2 WATER MOLECULES
       qxtte_a(:,idt_lg_H2O) = qxtte_a(:,idt_lg_H2O) - 2._dp*CH4_te_lg(:)
! qqq-
    END IF

    IF (i_H2O_feedback == 2) THEN
       CALL mtend_add_g(my_handle, mtend_id_q, px=q_te)
    END IF
#endif

    ! CLASS++
    IF (l_lg_class) &
         CALL class_integrate_lg( &
         temp_lg(:), press_3d_lg(:), spechum_lg(:) )
    ! CLASS--

    ! ISO++
    IF (l_lg_iso_C .OR. l_lg_iso_H) &
         CALL iso_integrate_lg( &
         temp_lg(:), press_3d_lg(:), spechum_lg(:), CH4_te_lg(:) )
    ! ISO--

    ! CLEAN MEMORY
    DEALLOCATE(temp)       ; NULLIFY(temp)
    DEALLOCATE(spechum)    ; NULLIFY(spechum)
    DEALLOCATE(temp_lg)    ; NULLIFY(temp_lg)
    DEALLOCATE(spechum_lg) ; NULLIFY(spechum_lg)
    DEALLOCATE(CH4_lg)     ; NULLIFY(CH4_lg)
    DEALLOCATE(CH4_te_lg)  ; NULLIFY(CH4_te_lg)
    !
    IF (i_H2O_feedback == 2) THEN
       DEALLOCATE(q_te)    ; NULLIFY(q_te)
       DEALLOCATE(q_te_lg) ; NULLIFY(q_te_lg)
    END IF

#endif
!!#D attila -

  END SUBROUTINE ch4_global_end
  ! ====================================================================

!!#D attila +
#ifdef ECHAM5
  ! === PRIVATE ========================================================
  SUBROUTINE class_integrate_lg(temp, press, spechum)

    USE messy_main_timer,         ONLY: time_step_len
    USE messy_main_tracer_mem_bi, ONLY: qxtte_a, qxtm1_a, NCELL
    USE messy_main_tools,         ONLY: PTR_1D_ARRAY

    IMPLICIT NONE

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(IN) :: temp
    REAL(DP), DIMENSION(:), INTENT(IN) :: press
    REAL(DP), DIMENSION(:), INTENT(IN) :: spechum

    ! LOCAL
    INTEGER                                     :: je, ja 
    INTEGER                                     :: idt
    TYPE(PTR_1D_ARRAY), DIMENSION(:,:), POINTER :: CH4c    => NULL()
    TYPE(PTR_1D_ARRAY), DIMENSION(:,:), POINTER :: CH4c_te => NULL()

    ! allocate space for tracers and tendencies
    ! (need to be stored individually for adjustment)
    ALLOCATE(CH4c(i_lg_nclass_emis_age(1), i_lg_nclass_emis_age(2)))
    ALLOCATE(CH4c_te(i_lg_nclass_emis_age(1), i_lg_nclass_emis_age(2)))

    emis_lg: DO je=1, i_lg_nclass_emis_age(1)
       age_lg: DO ja=1,  i_lg_nclass_emis_age(2)

          ! allocate space for tracers and tendencies
          ! (need to be stored individually for adjustment)
          ALLOCATE(CH4c(je,ja)%ptr(NCELL))
          CH4c(je,ja)%ptr(:) = 0.0_dp
          ALLOCATE(CH4c_te(je,ja)%ptr(NCELL))
          CH4c_te(je,ja)%ptr(:) = 0.0_dp

          ! calculate individual start values
!qqq+ NOTE: TENDENCY NOT YET IMPLEMENTED FOR LG TRACERS ...
!#ifndef MESSYTENDENCY
          idt = idt_lg_class(je,ja)
          CH4c(je,ja)%ptr(:) = qxtm1_a(:,idt) &
               + qxtte_a(:,idt) * time_step_len
!#else
!         CALL mtend_get_start_l(mtend_id_tracer, v0=CH4c_lg(je,ja)%ptr &
!              , idt=idt_lg_class(je,ja))
!#endif
!qqq-

          ! integrate chemistry and store individual tendencies
          CALL ch4_integrate(        &
               CH4c_te(je,ja)%ptr(:) &
               , CH4c(je,ja)%ptr(:)  &
               , ptr_lg_OH(:)        &
               , ptr_lg_O1D(:)       &
               , ptr_lg_Cl(:)        &
               , ptr_lg_jCH4(:)      &
               , temp(:)             &
               , press(:)            &
               , spechum(:)          &
               )

       END DO age_lg
    END DO emis_lg

    CALL class_age_move_lg(CH4c, CH4c_te)

    IF (l_lg_adj_tend) CALL class_adj_tend_lg(CH4c, CH4c_te)

    ! ADD TRACER TENDENCIES
    DO je=1, i_lg_nclass_emis_age(1)
       DO ja=1,  i_lg_nclass_emis_age(2)
!qqq+ NOTE: TENDENCY NOT YET IMPLEMENTED FOR LG TRACER 
!#ifndef MESSYTENDENCY
          idt = idt_lg_class(je,ja)
          qxtte_a(:,idt) = &
               qxtte_a(:,idt) + CH4c_te(je,ja)%ptr(:)
!#else
!         CALL mtend_add_l(my_handle, mtend_id_tracer &
!             , px=CH4c_te_lg(je,ja)%ptr, idt=idt_lg_class(je,ja))
!#endif
          
          ! not longer needed, clean memory
          DEALLOCATE(CH4c(je,ja)%ptr)   ; NULLIFY(CH4c(je,ja)%ptr)
          DEALLOCATE(CH4c_te(je,ja)%ptr); NULLIFY(CH4c_te(je,ja)%ptr)
       END DO
    END DO

    ! CLEAN MEMORY
    DEALLOCATE(CH4c)      ; NULLIFY(CH4c)
    DEALLOCATE(CH4c_te)   ; NULLIFY(CH4c_te)

  END SUBROUTINE class_integrate_lg
  ! ====================================================================

  ! === PRIAVE =========================================================
  SUBROUTINE class_age_move_lg(CH4c, CH4c_te)

    ! UPDATE AGE CLASSES

    USE messy_main_tools,         ONLY: PTR_1D_ARRAY
    USE messy_main_timer,         ONLY: time_step_len

    USE messy_main_tracer_mem_bi, ONLY: qxt_a, qxtm1_a, qxtf_a, qxtte_a
    USE messy_main_data_bi,       ONLY: eps

    IMPLICIT NONE

    ! I/O
    TYPE(PTR_1D_ARRAY), DIMENSION(:,:), POINTER :: CH4c
    TYPE(PTR_1D_ARRAY), DIMENSION(:,:), POINTER :: CH4c_te

    ! LOCAL
    INTEGER  :: ja, je ! age and emission class counter
    REAL(DP) :: sro    ! scale for 'remove own'
    INTEGER  :: idt, idt2

    IF (.NOT. (l_lg_mv_agecl_now .AND. (i_lg_nclass_emis_age(2) > 1))) RETURN

    SELECT CASE(i_lg_ageing)

    CASE(0,1)

       DO ja = i_lg_nclass_emis_age(2), 2, -1
          IF (ja == i_lg_nclass_emis_age(2)) THEN
             sro = 0.0_dp
          ELSE
             sro = -1.0_dp
          END IF
          DO je=1, i_lg_nclass_emis_age(1)
             ! tendency to remove entire own class  (not for oldest!)
             ! + tendency to add entire younger class
             CH4c_te(je,ja)%ptr(:) = &
                  f_lg(ja) * ( &
                  sro * ( ( CH4c(je,ja)%ptr(:) &
                  + CH4c_te(je,ja)%ptr(:) * time_step_len )/time_step_len ) ) &
                  + f_lg(ja-1) * ( ( CH4c(je,ja-1)%ptr(:) &
                  + CH4c_te(je,ja-1)%ptr(:) * time_step_len )/time_step_len )
          END DO
       END DO
       ja = 1
       DO je=1, i_lg_nclass_emis_age(1)
          CH4c_te(je,ja)%ptr(:) = &
               f_lg(ja) * ( &
               -1.0_dp * ( ( CH4c(je,ja)%ptr(:) &
               + CH4c_te(je,ja)%ptr(:) * time_step_len )/time_step_len ) &
               )
       END DO

    CASE(2)
       ! RESET TRACERS COMPLETELY
       DO ja = i_lg_nclass_emis_age(2), 2, -1
          IF (ja == i_lg_nclass_emis_age(2)) THEN
             sro = 0.0_dp
          ELSE
             sro = -1.0_dp
          END IF
          DO je=1, i_lg_nclass_emis_age(1)
             idt  = idt_lg_class(je,ja)
             idt2 = idt_lg_class(je,ja-1)
             IF (ASSOCIATED(qxt_a)) &
                  qxt_a(:,idt) = qxt_a(:,idt) + &
                  sro * qxt_a(:,idt) + &
                  qxt_a(:,idt2)
             IF (ASSOCIATED(qxtm1_a)) &
                  qxtm1_a(:,idt) = (1.0_dp - eps) *  &
                  qxt_a(:,idt)
             IF (ASSOCIATED(qxtf_a)) &
                  qxtf_a(:,idt) = qxt_a(:,idt) + &
                  eps*(qxtm1_a(:,idt) - 2.0_dp*qxt_a(:,idt))
             IF (ASSOCIATED(qxtte_a)) &
                  qxtte_a(:,idt) = qxtte_a(:,idt) + & 
                  sro * qxtte_a(:,idt) + &
                  qxtte_a(:,idt2)
          END DO
       END DO
       ja = 1
       DO je=1, i_lg_nclass_emis_age(1)
          idt  = idt_lg_class(je,ja)
          IF (ASSOCIATED(qxt_a))   qxt_a(:,idt)   = 0.0_dp
          IF (ASSOCIATED(qxtm1_a)) qxtm1_a(:,idt) = 0.0_dp
          IF (ASSOCIATED(qxtf_a))  qxtf_a(:,idt)  = 0.0_dp
          IF (ASSOCIATED(qxtte_a)) qxtte_a(:,idt) = 0.0_dp
       END DO

    END SELECT

  END SUBROUTINE class_age_move_lg
  ! =========================================================================

  ! === PRIVATE =============================================================
  SUBROUTINE class_adj_tend_lg(CH4c, CH4c_te)

    USE messy_main_timer,         ONLY: time_step_len
    USE messy_main_tracer_mem_bi, ONLY: qxtte_a, qxtm1_a &
                                      , NCELL
    USE messy_main_tools,         ONLY: PTR_1D_ARRAY

    IMPLICIT NONE

    ! I/O
    TYPE(PTR_1D_ARRAY), DIMENSION(:,:), POINTER :: CH4c
    TYPE(PTR_1D_ARRAY), DIMENSION(:,:), POINTER :: CH4c_te
    
    ! LOCAL
    INTEGER                         :: je, ja
    INTEGER                         :: idt
    REAL(DP), DIMENSION(:), POINTER :: CH4_m1_sum   => NULL() 
    REAL(DP), DIMENSION(:), POINTER :: CH4_te_sum   => NULL() 
    REAL(DP), DIMENSION(:), POINTER :: tmp => NULL()

    ! --- ADJUST TENDENCIES TO MASTER TRACER .......................
    ! --- 1) MEMORY
    ALLOCATE(CH4_m1_sum(NCELL))
    CH4_m1_sum(:) = 0.0_dp
    ALLOCATE(CH4_te_sum(NCELL))
    CH4_te_sum(:) = 0.0_dp
    ALLOCATE(tmp(NCELL))
    tmp(:) = 0.0_dp
    !
    sadj_lg(:) = 0.0_dp
    !
    ! --- 2) calculate SUM
    DO je=1, i_lg_nclass_emis_age(1)
       DO ja=1,  i_lg_nclass_emis_age(2)
          idt = idt_lg_class(je,ja)
          CH4_m1_sum(:) = CH4_m1_sum(:) + qxtm1_a(:,idt)
          CH4_te_sum(:) = CH4_te_sum(:) + qxtte_a(:,idt) + &
               CH4c_te(je,ja)%ptr(:)   ! not yet added !!!
       END DO
    END DO
    !
    ! --- 3) calculate scaling (weight) factor
    ! Note: 'CH4' tendency has already been added above to master tracer ...
    idt = idt_lg_CH4 ! master tracer
    CALL sca_tend(qxtm1_a(:,idt), qxtte_a(:,idt) &
         , CH4_m1_sum(:), CH4_te_sum(:), time_step_len, sadj_lg(:))
    !
    ! --- 4) apply scaling / add additional tendency for adjustment
    DO je=1, i_lg_nclass_emis_age(1)
       DO ja=1,  i_lg_nclass_emis_age(2)
          idt = idt_lg_class(je,ja)
          CALL adj_tend(qxtm1_a(:,idt)  &
               , qxtte_a(:,idt)         &
               , sadj_lg(:)             &
               , time_step_len          &
               , tmp(:)                 &
               )
          CH4c_te(je,ja)%ptr(:) = CH4c_te(je,ja)%ptr(:) + tmp(:)
       END DO
    END DO
    ! --- END ADJUST TENDENCIES TO MASTER TRACER ....................

    ! CLEAN MEMORY
    DEALLOCATE(tmp)       ; NULLIFY(tmp)
    DEALLOCATE(CH4_m1_sum); NULLIFY(CH4_m1_sum)
    DEALLOCATE(CH4_te_sum); NULLIFY(CH4_te_sum)

  END SUBROUTINE class_adj_tend_lg
  ! ====================================================================
#endif
!!#D attila -

!ISO++
  ! ####################################################################
  ! PRIVATE SUBROUTINES for isotopologue tracers
  ! ####################################################################

  ! === PRIVATE ========================================================
  SUBROUTINE iso_integrate_gp(temp, press, spechum, CH4_te) 

    USE messy_main_constants_mem,  ONLY: M_air, M_H2O
    USE messy_main_timer,          ONLY: time_step_len
    USE messy_main_grid_def_mem_bi,ONLY: kproma, nlev, jrow 
    USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte, pxtm1 => qxtm1
    USE messy_main_tools,          ONLY: PTR_2D_ARRAY

    IMPLICIT NONE

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: temp
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: press
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: spechum
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: CH4_te


    ! LOCAL
    INTEGER                           :: idt
    INTEGER                           :: ii 

    TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER :: CH4c    => NULL()
    TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER :: CH4c_te => NULL()

    ! HDO++
    REAL(DP),               PARAMETER :: scvmr           = M_air/M_H2O
    REAL(DP), DIMENSION(:,:), POINTER :: HDO_te => NULL() 

    ! HDO--

    ! Allocate
    ALLOCATE(CH4c(size(iso_id_gp)))
    ALLOCATE(CH4c_te(size(iso_id_gp)))
    IF (HDO_ex_gp) ALLOCATE(HDO_te(kproma, nlev))

    ! For every subtracer
    DO ii=1,size(iso_id_gp)

       ALLOCATE(CH4c(ii)%ptr(kproma, nlev))
       ALLOCATE(CH4c_te(ii)%ptr(kproma, nlev))

    ! Start values
#ifndef MESSYTENDENCY
       idt = idt_gp_iso(ii)
       
       CH4c(ii)%ptr(:,:) = pxtm1(_RI_X_ZN_(1:kproma,:,idt)) &       
            + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
#else
       CALL mtend_get_start_l(idt_gp_iso(ii), v0=CH4c(ii)%ptr(:,:)) 
       
#endif
       ! integrate chemistry and store individual tendencies
       CH4c_te(ii)%ptr(:,:) = 0.0_dp

       ! ch4_integrate for isotopolouge iso_id_gp(ii)
       CALL ch4_integrate(&
            CH4c_te(ii)%ptr(1:kproma,:)            &
            , CH4c(ii)%ptr(1:kproma,:)             &
            , ptr_gp_OH(_RI_XYZ__(1:kproma,jrow,:))    &
            , ptr_gp_O1D(_RI_XYZ__(1:kproma,jrow,:))   &
            , ptr_gp_Cl(_RI_XYZ__(1:kproma,jrow,:))    &
            , ptr_gp_jCH4(_RI_XYZ__(1:kproma,jrow,:))  &
            , temp(:,:)            &
            , press(:,:)  &
            , spechum(:,:)         &
            , iso_id=iso_id_gp(ii)    &
            )
    END DO
          
    ! HDO++
    IF (HDO_ex_gp) THEN
       IF ( l_ef_re ) THEN
          ! Additionally we have to add the correction for the D, that stays
          ! in HD reservoirs, see Eichinger et al. (2015) Part 1
          !    hdo_te(1:kproma,:) = -1._dp * ( (2._dp*CH3D_te(1:kproma,:)) - &
          !     (6.32e-5_dp*CH4_te(1:kproma,:)) ) &
          !          / (scvmr_D * (1.0_dp / (1.0_dp - HDO_F(1:kproma,:)))**2)
          ! In current syntax:
          HDO_te(1:kproma,:) = &
               -1._dp * ( (CH4c_te(i_gp_CH4_D1)%ptr(1:kproma,:)) &
               - (6.32e-5_dp*CH4_te(1:kproma,:)) )
       ELSE
       ! One CH3D molecule reacts to one HDO molecule
       ! => calculate tendency of HDO
       ! => add tendency to tracer (see below)
          HDO_te(1:kproma,:) = &
               (-1._dp*CH4c_te(i_gp_CH4_D1)%ptr(1:kproma,:)) 
       END IF

    END IF
    ! HDO--

    ! CHECK FOR CONSISTENCY OF SUBTRACERS (pairwise)
    IF (l_gp_adj_tend) THEN
       DO ii=1,size(iso_id_gp),2
          CALL iso_adj_tend_gp(CH4c(ii:ii+1), &
               CH4c_te(ii:ii+1), &
               idt_gp_iso(ii:ii+1) )
       END DO
    END IF
    
    ! ADD TRACER TENDENCIES
    DO ii=1,size(iso_id_gp)
#ifndef MESSYTENDENCY
       idt = idt_gp_iso(ii)
       pxtte(_RI_X_ZN_(1:kproma,:,idt)) = pxtte(_RI_X_ZN_(1:kproma,:,idt)) + &
            CH4c_te(ii)%ptr(1:kproma,:)
#else
       CALL mtend_add_l(my_handle, idt_gp_iso(ii), px=CH4c_te(ii)%ptr(:,:)) 
#endif

       ! clean memory
       DEALLOCATE(CH4c(ii)%ptr)   ; NULLIFY(CH4c(ii)%ptr)
       DEALLOCATE(CH4c_te(ii)%ptr); NULLIFY(CH4c_te(ii)%ptr)
    END DO

    ! ADD HDO
    IF (HDO_ex_gp) THEN
#ifndef MESSYTENDENCY
       idt = idt_gp_HDO
       pxtte(_RI_X_ZN_(1:kproma,:,idt)) = pxtte(_RI_X_ZN_(1:kproma,:,idt)) + HDO_te(1:kproma,:)
#else
       CALL mtend_add_l(my_handle, idt_gp_HDO, px=HDO_te(:,:))
#endif
    END IF

 
    ! CLEAN MEMORY
    DEALLOCATE(CH4c)      ; NULLIFY(CH4c)
    DEALLOCATE(CH4c_te)   ; NULLIFY(CH4c_te)

    IF (HDO_ex_gp) DEALLOCATE(HDO_te)  ; NULLIFY(HDO_te)

  END SUBROUTINE iso_integrate_gp
  ! =========================================================================

  ! === PRIVATE =============================================================
  SUBROUTINE iso_adj_tend_gp(CH4c, CH4c_te, idt_gp_iso_adj)

    USE messy_main_grid_def_mem_bi, ONLY: kproma, nlev, jrow
    USE messy_main_timer,           ONLY: time_step_len
    USE messy_main_tracer_mem_bi,   ONLY: pxtte => qxtte, pxtm1 => qxtm1 
    USE messy_main_tools,           ONLY: PTR_2D_ARRAY
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer,          ONLY: get_tracer
    USE messy_main_constants_mem,   ONLY: STRLEN_MEDIUM

    IMPLICIT NONE
 
   ! I/O
    TYPE(PTR_2D_ARRAY), DIMENSION(:)          :: CH4c    ! tracer concentrations
    TYPE(PTR_2D_ARRAY), DIMENSION(:)          :: CH4c_te ! tracer tendencies
    INTEGER, DIMENSION(:)                     :: idt_gp_iso_adj    ! list of tracer IDs

    ! LOCAL
    INTEGER                             :: idt
    INTEGER                             :: ii
    REAL(DP), DIMENSION(:,:),   POINTER :: CH4_m1_sum   => NULL() 
    REAL(DP), DIMENSION(:,:),   POINTER :: CH4_te_sum   => NULL() 
    REAL(DP), DIMENSION(:,:),   POINTER :: tmp => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: sadj_iso_gp => NULL()
    INTEGER                             :: status
    CHARACTER(LEN=STRLEN_MEDIUM)        :: sub

    ! --- ADJUST TENDENCIES TO MASTER TRACER .......................
    ! --- 1) MEMORY
    ALLOCATE(CH4_m1_sum(kproma, nlev))
    CH4_m1_sum(:,:) = 0.0_dp
    ALLOCATE(CH4_te_sum(kproma, nlev))
    CH4_te_sum(:,:) = 0.0_dp
    ALLOCATE(tmp(kproma, nlev))
    tmp(:,:) = 0.0_dp
    !
    ! get name of tracer by its ID
    CALL get_tracer(status, GPTRSTR, idt_gp_iso_adj(1) &
         , subname=sub)
    ! 
    ! Choose weight channel object according to subtracer
    SELECT CASE(sub)
      CASE('12C','13C') 
         sadj_iso_gp => sadj_iso_C_gp
      CASE('D0 ','D1 ')
         sadj_iso_gp => sadj_iso_H_gp
    END SELECT
    sadj_iso_gp(_RI_XYZ__(:,jrow,:)) = 0.0_dp

    !
    ! --- 2) calculate SUM
    DO ii=1,size(idt_gp_iso_adj)
       idt = idt_gp_iso_adj(ii)
       CH4_m1_sum(:,:) = CH4_m1_sum(:,:) + pxtm1(_RI_X_ZN_(1:kproma,:,idt))
       CH4_te_sum(:,:) = CH4_te_sum(:,:) + pxtte(_RI_X_ZN_(1:kproma,:,idt)) + &
               CH4c_te(ii)%ptr(:,:)   ! not yet added !!!
    END DO
    !
    ! --- 3) calculate scaling (weight) factor
    ! Note: 'CH4' tendency has already been added above to master tracer ...
    idt = idt_gp_CH4 ! master tracer
    CALL sca_tend(pxtm1(_RI_X_ZN_(1:kproma,:,idt))  &  ! master tracer
                 , pxtte(_RI_X_ZN_(1:kproma,:,idt)) &  ! master tendency
                 , CH4_m1_sum(:,:)         &  ! sum of subtracers
                 , CH4_te_sum(:,:)         &  ! sum of subtracer tendencies
                 , time_step_len           &  ! time step
                 , sadj_iso_gp(_RI_XYZ__(1:kproma,jrow,:)))! weights (OUTPUT)

    ! --- 4) apply scaling / add additional tendency for adjustment
    DO ii=1,2
       idt = idt_gp_iso_adj(ii)
       CALL adj_tend(pxtm1(_RI_X_ZN_(1:kproma,:,idt))  &
            , pxtte(_RI_X_ZN_(1:kproma,:,idt))         &
            , sadj_iso_gp(_RI_XYZ__(1:kproma,jrow,:))   &
            , time_step_len                   &
            , tmp(:,:)                        &
            )
       CH4c_te(ii)%ptr(:,:) = CH4c_te(ii)%ptr(:,:) + tmp(:,:)
    END DO

    ! --- END ADJUST TENDENCIES TO MASTER TRACER ....................

    ! CLEAN MEMORY
    DEALLOCATE(tmp)       ; NULLIFY(tmp)
    DEALLOCATE(CH4_m1_sum); NULLIFY(CH4_m1_sum)
    DEALLOCATE(CH4_te_sum); NULLIFY(CH4_te_sum)

  END SUBROUTINE iso_adj_tend_gp
  ! ====================================================================

!!#D attila +
#ifdef ECHAM5
  ! === PRIVATE ========================================================
  SUBROUTINE iso_integrate_lg(temp, press, spechum, CH4_te)

    USE messy_main_constants_mem, ONLY: M_air, M_H2O
    USE messy_main_timer,         ONLY: time_step_len
    USE messy_main_tracer_mem_bi, ONLY: qxtte_a, qxtm1_a &
                                         , NCELL
    USE messy_main_tools,         ONLY: PTR_1D_ARRAY

    IMPLICIT NONE

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(IN) :: temp
    REAL(DP), DIMENSION(:), INTENT(IN) :: press
    REAL(DP), DIMENSION(:), INTENT(IN) :: spechum
    REAL(DP), DIMENSION(:), INTENT(IN) :: CH4_te

    ! LOCAL
    INTEGER                                   :: idt
    INTEGER                                   :: ii

    TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER :: CH4c    => NULL()
    TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER :: CH4c_te => NULL()

    ! HDO++
    REAL(DP),               PARAMETER :: scvmr = M_air/M_H2O
    REAL(DP), DIMENSION(:), POINTER   :: HDO_te => NULL() 
    ! HDO--

    ! Allocate
    ALLOCATE(CH4c(size(iso_id_lg)))
    ALLOCATE(CH4c_te(size(iso_id_lg)))

    ! For every subtracer
    DO ii=1,size(iso_id_lg)
       
       ALLOCATE(CH4c(ii)%ptr(NCELL))
       ALLOCATE(CH4c_te(ii)%ptr(NCELL))
       IF (HDO_ex_lg) ALLOCATE(HDO_te(NCELL))

       ! Start values
!qqq+ NOTE: TENDENCY NOT YET IMPLEMENTED FOR LG TRACERS ...
!#ifndef MESSYTENDENCY
       idt = idt_lg_iso(ii)
       CH4c(ii)%ptr(:) = qxtm1_a(:,idt) &
            + qxtte_a(:,idt) * time_step_len
!#else
!         CALL mtend_get_start_l(mtend_id_tracer, v0=CH4c(ii)%ptr(:) &
!              , idt=idt_lg_iso(ii))
!       
!#endif

       ! integrate chemistry and store individual tendencies
       CH4c_te(ii)%ptr(:) = 0.0_dp

       CALL ch4_integrate(        &
            CH4c_te(ii)%ptr(:) &
            , CH4c(ii)%ptr(:)  &
            , ptr_lg_OH(:)        &
            , ptr_lg_O1D(:)       &
            , ptr_lg_Cl(:)        &
            , ptr_lg_jCH4(:)      &
            , temp(:)             &
            , press(:)            &
            , spechum(:)          &
            , iso_id=iso_id_lg(ii) &
            )
    END DO
          
    ! HDO++
    IF (HDO_ex_lg) THEN

       IF ( l_ef_re ) THEN
          ! Additionally we have to add the correction for the D, that stays
          ! in HD reservoirs, see Eichinger et al. (2015) Part 1
          !    hdo_te(1:kproma,:) = -1._dp * ( (2._dp*CH3D_te(1:kproma,:)) - &
          !     (6.32e-5_dp*CH4_te(1:kproma,:)) ) &
          !          / (scvmr_D * (1.0_dp / (1.0_dp - HDO_F(1:kproma,:)))**2)
          HDO_te(:) = &
               -1._dp * ( (CH4c_te(i_lg_CH4_D1)%ptr(:)) &
               - (6.32e-5_dp*CH4_te(:)) )
       ELSE
       ! One CH3D molecule reacts to one HDO molecule
       ! => calculate tendency of HDO
       ! => add tendency to tracer (see below)
          HDO_te(:) = (-1._dp*CH4c_te(i_lg_CH4_D1)%ptr(:))
       END IF
    END IF
    ! HDO--

    ! CHECK FOR CONSISTENCY OF SUBTRACERS (pairwise)
    IF (l_lg_adj_tend) THEN
       DO ii=1,size(iso_id_lg),2
          CALL iso_adj_tend_lg(CH4c(ii:ii+1), &
               CH4c_te(ii:ii+1), &
               idt_lg_iso(ii:ii+1) )
       END DO
    END IF

    ! ADD TRACER TENDENCIES
    !qqq+ NOTE: TENDENCY NOT YET IMPLEMENTED FOR LG TRACER 
    DO ii=1,size(iso_id_lg)
!#ifndef MESSYTENDENCY
       idt = idt_lg_iso(ii)
       qxtte_a(:,idt) = &
            qxtte_a(:,idt) + CH4c_te(ii)%ptr(:)
!#else  
!        CALL mtend_add_l(my_handle, mtend_id_tracer, px=CH4c_te(ii)%ptr(:), &
!             idt=idt_lg_iso(ii))
!#endif  
       DEALLOCATE(CH4c(ii)%ptr)   ; NULLIFY(CH4c(ii)%ptr)
       DEALLOCATE(CH4c_te(ii)%ptr); NULLIFY(CH4c_te(ii)%ptr)
  
    END DO

    ! ADD HDO
    IF (HDO_ex_lg) THEN
!#ifndef MESSYTENDENCY
       idt = idt_lg_HDO
       qxtte_a(:,idt) = qxtte_a(:,idt) + HDO_te(:)
! #else
!        CALL mtend_add_l(my_handle, mtend_id_tracer, &
!            px=HDO_te(:))
! #endif
    END IF
    
    ! CLEAN MEMORY
    DEALLOCATE(CH4c)      ; NULLIFY(CH4c)
    DEALLOCATE(CH4c_te)   ; NULLIFY(CH4c_te)

    IF (HDO_ex_lg) DEALLOCATE(HDO_te)  ; NULLIFY(HDO_te)

  END SUBROUTINE iso_integrate_lg
  ! =========================================================================

  ! === PRIVATE =============================================================
  SUBROUTINE iso_adj_tend_lg(CH4c, CH4c_te, idt_lg_iso_adj)

    USE messy_main_timer,         ONLY: time_step_len
    USE messy_main_tracer_mem_bi, ONLY: qxtte_a, qxtm1_a &
                                      , NCELL, LGTRSTR
    USE messy_main_tools,         ONLY: PTR_1D_ARRAY
    USE messy_main_tracer,        ONLY: get_tracer
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM

    IMPLICIT NONE

    ! I/O
    TYPE(PTR_1D_ARRAY), DIMENSION(:)  :: CH4c
    TYPE(PTR_1D_ARRAY), DIMENSION(:)  :: CH4c_te
    INTEGER, DIMENSION(:)             :: idt_lg_iso_adj 

    ! LOCAL
    INTEGER                             :: idt
    INTEGER                             :: ii
    INTEGER                             :: status
    REAL(DP), DIMENSION(:), POINTER     :: CH4_m1_sum   => NULL() 
    REAL(DP), DIMENSION(:), POINTER     :: CH4_te_sum   => NULL() 
    REAL(DP), DIMENSION(:), POINTER     :: tmp => NULL()
    REAL(DP), DIMENSION(:), POINTER     :: sadj_iso_lg => NULL()
    CHARACTER(LEN=STRLEN_MEDIUM)        :: sub

    ! --- ADJUST TENDENCIES TO MASTER TRACER .......................
    ! --- 1) MEMORY
    ALLOCATE(CH4_m1_sum(NCELL))
    CH4_m1_sum(:) = 0.0_dp
    ALLOCATE(CH4_te_sum(NCELL))
    CH4_te_sum(:) = 0.0_dp
    ALLOCATE(tmp(NCELL))
    tmp(:) = 0.0_dp
    !
    ! get name of tracer by its ID
    CALL get_tracer(status, LGTRSTR, idt_lg_iso_adj(1) &
         , subname=sub)
    ! 
    SELECT CASE(sub)
      CASE('12C','13C') 
         sadj_iso_lg => sadj_iso_C_lg
      CASE('D0 ','D1 ')
         sadj_iso_lg => sadj_iso_H_lg
    END SELECT
    sadj_iso_lg(:) = 0.0_dp
    !
    ! --- 2) calculate SUM
    DO ii=1,size(idt_lg_iso_adj)
       idt = idt_lg_iso_adj(ii)
       CH4_m1_sum(:) = CH4_m1_sum(:) + qxtm1_a(:,idt)
       CH4_te_sum(:) = CH4_te_sum(:) + qxtte_a(:,idt) + &
            CH4c_te(ii)%ptr(:)   ! not yet added !!!
    END DO
    !
    ! --- 3) calculate scaling (weight) factor
    ! Note: 'CH4' tendency has already been added above to master tracer ...
    idt = idt_lg_CH4 ! master tracer
    CALL sca_tend(qxtm1_a(:,idt), qxtte_a(:,idt) &
         , CH4_m1_sum(:), CH4_te_sum(:), time_step_len &
         , sadj_iso_lg(:))
    !
    ! --- 4) apply scaling / add additional tendency for adjustment
    DO ii=1,2
       idt = idt_lg_iso_adj(ii)
       CALL adj_tend(qxtm1_a(:,idt)  &
            , qxtte_a(:,idt)         &
            , sadj_iso_lg(:)         &
            , time_step_len          &
            , tmp(:)                 &
            )
       CH4c_te(ii)%ptr(:) = CH4c_te(ii)%ptr(:) + tmp(:)
    END DO
    ! --- END ADJUST TENDENCIES TO MASTER TRACER ....................

    ! CLEAN MEMORY
    DEALLOCATE(tmp)       ; NULLIFY(tmp)
    DEALLOCATE(CH4_m1_sum); NULLIFY(CH4_m1_sum)
    DEALLOCATE(CH4_te_sum); NULLIFY(CH4_te_sum)

  END SUBROUTINE iso_adj_tend_lg
  ! ====================================================================
#endif
!!#D attila -
!ISO--

  ! ====================================================================
  SUBROUTINE ch4_free_memory

    IF (ASSOCIATED(idt_gp_class)) THEN
       DEALLOCATE(idt_gp_class)
       NULLIFY(idt_gp_class)
    END IF

    IF (ASSOCIATED(f_gp)) THEN
       DEALLOCATE(f_gp)
       NULLIFY(f_gp)
    END IF

    !ISO++
    IF (ASSOCIATED(idt_gp_iso)) THEN
       DEALLOCATE(idt_gp_iso)
       NULLIFY(idt_gp_iso)
    END IF

    IF (ASSOCIATED(iso_id_gp)) THEN
       DEALLOCATE(iso_id_gp)
       NULLIFY(iso_id_gp)
    END IF
    !ISO--

    IF (ASSOCIATED(idt_lg_class)) THEN
       DEALLOCATE(idt_lg_class)
       NULLIFY(idt_lg_class)
    END IF

    IF (ASSOCIATED(f_lg)) THEN
       DEALLOCATE(f_lg)
       NULLIFY(f_lg)
    END IF

    ! ISO++
    IF (ASSOCIATED(idt_lg_iso)) THEN
       DEALLOCATE(idt_lg_iso)
       NULLIFY(idt_lg_iso)
    END IF

    IF (ASSOCIATED(iso_id_lg)) THEN
       DEALLOCATE(iso_id_lg)
       NULLIFY(iso_id_lg)
    END IF
    ! ISO--

  END SUBROUTINE ch4_free_memory
  ! ====================================================================

  ! ####################################################################
  ! PRIVATE SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE ch4_read_nml_cpl(status, iou)
   
    ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close
    ! MESSy/BMIL
#ifdef ECHAM5
    USE messy_main_tracer_mem_bi, ONLY: NGCELL
#endif

    IMPLICIT NONE
    
    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ i_H2O_feedback &
         , L_GP, c_gp_OH, c_gp_O1D, c_gp_Cl, c_gp_jCH4 &
         , L_LG, c_lg_OH, c_lg_O1D, c_lg_Cl, c_lg_jCH4 &
         ! CLASS++
         , i_gp_nclass_emis_age, i_lg_nclass_emis_age  &
         , l_gp_adj_tend, l_lg_adj_tend                &
         , i_gp_ageing, i_lg_ageing                    &
         , r_gp_age_cll, r_lg_age_cll                  &
         ! CLASS--
         ! ISO++
         , l_gp_iso_C, l_gp_iso_H                      &
         , l_lg_iso_C, l_lg_iso_H                      &
         ! ISO--
         , l_ef_re

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='ch4_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    ! SET DEFAULT VALUES
    r_gp_age_cll(:) = 30.44_dp  ! average month length in days
    r_lg_age_cll(:) = 30.44_dp  ! average month length in days

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! ### ADD HERE DIAGNOSTIC OUTPUT FOR LOG-FILE
    IF ((.NOT.L_GP).AND.(.NOT.L_LG)) THEN
       CALL warning_bi('GRIDPOINT AND LAGRANGIAN INTEGRATION'//&
            &' BOTH SWITCHED OFF!', substr)
    END IF

!!#D attila +
#ifdef ECHAM5
    IF (L_LG) THEN
       IF (NGCELL <= 0) THEN
          CALL warning_bi('LAGRANGIAN INTEGRATION REQUIRES ATTILA!'//&
               &'-> SWITCHING OFF LAGRANGIAN INTEGRATION FOR CH4', substr)
          L_LG = .FALSE.
       END IF
    END IF
#endif
!!#D attila -

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE ch4_read_nml_cpl
  ! ====================================================================

! **********************************************************************
END MODULE messy_ch4_si
! **********************************************************************
