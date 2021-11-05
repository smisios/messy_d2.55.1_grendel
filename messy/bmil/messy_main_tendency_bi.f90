#include "messy_main_ppd_bi.inc"

module messy_main_tendency_bi

#if defined(MESSYTENDENCY)

  ! This module performs three main tasks:
  !
  !    1. Compute initial value and update tendencies
  !
  !    2. Collect individual tendencies and write them into output channel
  !
  !    3. Provide interface for requesting process-tendencies in other submodels
  !
  !   Authors:
  !   Klaus Ketelsen, Month Year:
  !     - first implementation for MESSy1
  !
  !   Roland Eichinger and Patrick Joeckel, March 2012:
  !     - completely rewritten for MESSy2
  !
  !   Astrid Kerkweg
  !     - completely revised for usage with other basemodel than ECHAM (2018)
  !     - completely revised for usage in ICON (2019-2020)
  !     - tendency memory handling revised (2020)

  ! MESSy/SMCL
  USE messy_main_constants_mem,    ONLY: dp, INT_UNDEF, iouerr
  USE messy_main_tools,            ONLY: PTR_3D_ARRAY, PTR_1D_ARRAY_INT
  USE messy_main_channel,          ONLY: t_chaobj_cpl
  USE messy_main_timer,            ONLY: time_step_len
  USE messy_main_tendency
  ! MESSy/BMIL
  USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi &
                                       , error_bi
  USE messy_main_mpi_bi,           ONLY: p_pe, p_io
  USE messy_main_tracer_mem_bi,    ONLY: xtte, xtm1
  USE messy_main_grid_def_mem_bi,  ONLY: nlev
#ifdef ECHAM5
!WWWWW
!  USE messy_main_grid_def_bi,          ONLY: sqcst_2d
  USE mo_geoloc,                   ONLY: sqcst_2d
#endif
  USE messy_main_channel_mem,      ONLY: n_dom,dom_curid
#ifdef ICON
  ! needs to be parameter
  USE messy_main_bmluse_bi,        ONLY: max_dom
#endif

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTRINSIC :: NULL

#ifndef ICON
  INTEGER, PARAMETER :: max_dom = 1 ! needs to be parameter
#endif

  ! Number of prognostic variables (without tracers)
  INTEGER, PARAMETER :: NProgVmax = -8
  ! prognostic variable identifier in this generic submodel
  ! (must be consistent with pvarname, pvarunit, etc.!)
  INTEGER, PARAMETER, PUBLIC :: mtend_id_w  = -8
  INTEGER, PARAMETER, PUBLIC :: mtend_id_pp = -7
  INTEGER, PARAMETER, PUBLIC :: mtend_id_t  = -6
  INTEGER,            PUBLIC :: mtend_id_q  = -5
  INTEGER,            PUBLIC :: mtend_id_xl = -4
  INTEGER,            PUBLIC :: mtend_id_xi = -3
  INTEGER, PARAMETER, PUBLIC :: mtend_id_u  = -2
  INTEGER, PARAMETER, PUBLIC :: mtend_id_v  = -1

  CHARACTER(LEN=4), DIMENSION(NProgVmax:-1) :: pvarname = (/ &
       'w   ','pp  ','t   ','q   ','xl  ','xi  ','u   ','v   ' /)
  CHARACTER(LEN=7), DIMENSION(NProgVmax:-1) :: pvarunit = (/ &
       'm/s/s  ','Pa/s   ','K/s    ','kg/kg/s','kg/kg/s','kg/kg/s', &
       'm/s/s  ','m/s/s  ' /)

  INTEGER, PARAMETER, PUBLIC :: mtend_id_tracer  = 0

  ! BASEMODEL DEPENDENT LIST OF PROGNOSTIC VARIABLES
  TYPE T_PROGVARS
     TYPE(t_chaobj_cpl)                  :: var
     TYPE(t_chaobj_cpl)                  :: varte
     LOGICAL                             :: lexists  =  .FALSE.
     REAL(dp), DIMENSION(:,:,:), POINTER :: varm1ptr => NULL()
     REAL(dp), DIMENSION(:,:,:), POINTER :: varteptr => NULL()
     INTEGER                             :: reprid
  END TYPE T_PROGVARS

  TYPE(T_PROGVARS), DIMENSION(:), POINTER ::progvar => NULL()

  ! NAMELIST PARAMETERS
  LOGICAL,DIMENSION(max_dom) :: l_full_diag = .FALSE.
  LOGICAL,DIMENSION(max_dom) :: l_closure   = .FALSE.
  LOGICAL,DIMENSION(max_dom) :: l_clos_diag = .FALSE.

  !  MEMORY / POINTER REQUEST SOURCES
  INTEGER, PARAMETER :: IR_NML = 1  ! NAMELIST DIAG
  INTEGER, PARAMETER :: IR_FUL = 2  ! FULL DIAG
  INTEGER, PARAMETER :: IR_EXC = 3  ! EXCHANGE OF TENDENCIES BETWEEN SMs !
  !                                 ! (INTERNAL USE ONLY)
  INTEGER, PARAMETER :: IR_CLS = 4  ! CLOSURE DIAG (I_HANDLE_-SUM/DIFF)
  INTEGER, PARAMETER :: NREQ   = 4  ! NUMBER OF possible sources of tendency requests

  ! CORRESPONDING CHANNEL NAME POST-FIXES
  CHARACTER(LEN=*), DIMENSION(NREQ), PARAMETER :: chpf = (/ &
       '_diag',  &  ! namelist based diagnostic
       '_full',  &  ! full diagnostic output
       '_exch',  &  ! exchange of tend. between submodels (internal use only)
       '_clsr' /)   ! closure relevant diagnostics (I_HANDLE_-SUM/DIFF)

  ! NUMBER OF DIFFERENT CASES OF TENDENCY STORAGE
  ! 1 -> GENEREAL SINGLE TENDENCY STORAGE
  ! 2 -> NAMELIST REQUESTED: possible sum of processes
  INTEGER, PARAMETER :: IC_GEN = 1
  INTEGER, PARAMETER :: IC_NML = 2
  INTEGER, PARAMETER :: NCASES = 2

  INTEGER                     :: iou    ! I/O unit
  INTEGER                     :: status

  TYPE T_POINTERS
     !
     ! PROGNOSTIC VARIABLES
     ! ... to avoid IF(ASSOCIATED(...))
     LOGICAL, DIMENSION(:), POINTER            :: lvar   => NULL()
     TYPE(PTR_3D_ARRAY), DIMENSION(:),POINTER  :: pvar   => NULL()
     ! The following logical indicates, that the memory of this tendency
     ! is managed OUTSIDE of TENDENCY, i.e. values to the fields are assigned
     ! outside of TENDENCY and mtend_add is only called to add the tendency
     ! to the closure (SUM) field. This also includes the reset of the tendency
     ! to zero.
     ! This exemption is required for ICON, where the TENDENCY memory is used
     ! for convection and turbulence in order to avoid unnecessary double
     ! allocation of the memory as ICON requires tracer tendency fields for
     ! the individual processes
     ! NOTE: The lextern switch is only used in the startup phase, mainly to
     ! check that not two different submodels require modification of the same
     ! tendency and to check if the memory is already allocated.
     ! To simplify the logic in the time loop, lvar is set .FALSE. for
     ! tendencies with lextern=.TRUE.
     LOGICAL, DIMENSION(:), POINTER            :: lextern => NULL()
     !
  END TYPE T_POINTERS

  ! INTRINSIC PROCESS HANDLES / PROG. VARIABLES -----------------------------
  ! record of process / prognostic variable (incl. tracers) combination
  TYPE T_record
    ! FLAG, WHICH VARIABLES ARE MODIFIED BY THIS PROCESS
     LOGICAL, DIMENSION(:), POINTER      :: lvmod => NULL()
     ! ### SPECIAL TREATMENT FOR (GP) TRACERS
     ! DOES THIS PROCESS MODIFY ANY OF THE (GP) TRACERS?
     LOGICAL                             :: any_tracer = .FALSE.
     INTEGER                             :: ntrac_proc   = -1
     ! INDICES OF TRACERS MODIFIED BY THIS PROCESS
     INTEGER, DIMENSION(:), POINTER      :: idx => NULL()
     !
     ! MEMORY / POINTER MANAGEMENT
     !
     ! FOUR CASES:
     ! IC_GEN = 1  GENERAL SINGLE TENDENCY DIAG
     ! IC_NML = 2  NAMELIST DIAG
     TYPE(T_POINTERS), DIMENSION(NCASES) :: case
     !
  END TYPE T_record

  TYPE T_record_proc
     ! UNIQUE NAME OF PROCESS
     CHARACTER(LEN=32)                     :: name = '                                '
     ! UNIQUE ID OF PROCESS (dynamically determined)
     INTEGER                               :: handle = 0
     ! process per domain/patch
     TYPE(T_RECORD), DIMENSION(:), POINTER :: p
  END TYPE T_record_proc

  ! actual number of process handles
  INTEGER            :: N_HANDLE    = 0
  ! maximum number of handles
  INTEGER,PARAMETER  :: NMAX_HANDLE    = 100
  TYPE(T_record_proc), DIMENSION(NMAX_HANDLE) :: proc

  ! SPECIAL RECORDS FOR "summation of all tendencies" and
  ! calculation of difference to "total internal tendency"
  INTEGER :: I_HANDLE_SUM  = 0
  INTEGER :: I_HANDLE_DIFF = 0

  ! BaseModel Total Tendency
  ! for those submodels which do not use tendencies troughout the time step
  ! (ie COSMO and ICON) to check full change during time step vs. summed up
  ! tendencies.
  INTEGER :: I_HANDLE_BMTT  = 0

  ! ------------------------------------------------------------------------

  ! USER DEFINED (CPL-namelist) DIAGNOSTIC REQUESTS ------------------------
  TYPE T_TDIAG_IO
     CHARACTER(LEN=32)                    :: variable_name = ''
     CHARACTER(LEN=132)                   :: process_name = ''
  END TYPE T_TDIAG_IO

  TYPE T_TDIAG_IOD
     CHARACTER(LEN=32)                    :: variable_name = ''
     CHARACTER(LEN=32)                    :: domain_name   = ''
     INTEGER                              :: domain_idx    = INT_UNDEF
     CHARACTER(LEN=132)                   :: process_name  = ''
  END TYPE T_TDIAG_IOD

  TYPE T_TDIAG
     TYPE(T_TDIAG_IOD)                           :: io
     ! INDEX of prognostic variable in progvar
     INTEGER                                     :: progvaridx = INT_UNDEF
     INTEGER                                     :: ntclass
     TYPE(PTR_1D_ARRAY_INT),DIMENSION(:),POINTER :: iproc => NULL()
     INTEGER,DIMENSION(:),POINTER                :: nproc
     TYPE(PTR_3D_ARRAY),DIMENSION(:),POINTER     :: obj => NULL()
  END TYPE T_TDIAG

  INTEGER                     :: N_TDIAG
  INTEGER,PARAMETER           :: NMAX_TDIAG = 100

  TYPE(T_TDIAG_IO), DIMENSION(NMAX_TDIAG)       :: TDIAG  ! CPL
  TYPE(T_TDIAG), DIMENSION(:), POINTER          :: XTDIAG => NULL()
  INTEGER :: nump = 1

  ! Wind scaling
  ! default: horizontal wind vector is scaled with cos latitude
#ifdef ECHAM5
#define WINDSCALING 1
  logical            :: wind_is_scaled = .true.
#else
  logical            :: wind_is_scaled = .false.
#endif
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: lany_tracer

  ! ------------------------------------------------------------------------

  ! Interface section

  ! MAIN ENTRY POINTS (CALLED ONCE FROM BMIL)
  interface main_tendency_set_domain
     module procedure main_tendency_set_domain
  end interface

  interface main_tendency_initialize
     module procedure main_tendency_initialize
  end interface

  interface main_tendency_init_coupling
     module procedure main_tendency_init_coupling
  end interface

  interface main_tendency_global_end
     module procedure main_tendency_global_end
  end interface

  interface compute_eps_and_clear
     module procedure compute_eps_and_clear
  end interface

  interface main_tendency_free_memory
     module procedure main_tendency_free_memory
  end interface

  interface main_tendency_reset ! called from BML: scan1.f90)
     module procedure main_tendency_reset
  end interface

  ! ENTRY POINTS (CALLED FROM SMIL OF VARIOUS SUBMODELS)
  interface mtend_get_handle
     module procedure mtend_get_handle
  end interface

  interface mtend_register
     module procedure mtend_register
  end interface

  interface mtend_add_l
     module procedure mtend_add_l
  end interface

  interface mtend_add_g
     module procedure  mtend_add_g
  end interface

  interface mtend_get_start_l
     module procedure mtend_get_start_l
  end interface

  interface mtend_get_start_g
     module procedure  mtend_get_start_g
  end interface

  ! ENTRY POINT (CALLED FROM BML: physc.f90)
  interface mtend_set_sqcst_scal
     module procedure mtend_set_sqcst_scal
  end interface

  ! ENTRY POINT (CALLED FROM user in any feasible submodel)
  interface mtend_request
     module procedure mtend_request
  end interface

  ! op_ff_20160921+
  interface mtend_checkreg
     module procedure mtend_checkreg
  end interface
  ! op_ff_20160921-

  ! MAIN ENTRY POINTS (CALLED FROM MAIN_CONTROL)
  PUBLIC :: main_tendency_set_domain
  PUBLIC :: main_tendency_initialize
#ifdef COSMO
  PUBLIC :: main_tendency_new_tracer
#endif
  PUBLIC :: main_tendency_init_coupling
  PUBLIC :: main_tendency_global_end
  PUBLIC :: main_tendency_reset
  PUBLIC :: main_tendency_free_memory

  ! ENTRY POINTS CALLED FROM SMIL
  PUBLIC :: mtend_get_handle
  PUBLIC :: mtend_register
  PUBLIC :: mtend_get_start_l
  PUBLIC :: mtend_get_start_g
  PUBLIC :: mtend_add_l
  PUBLIC :: mtend_add_g

  PUBLIC :: mtend_request
  PUBLIC :: mtend_set_sqcst_scal

  ! op_ff_20160921+
  PUBLIC :: mtend_checkreg
  ! op_ff_20160921-

  !---------------------
  !  private subroutines
  !---------------------

  PRIVATE :: tendency_read_nml_cpl
  PRIVATE :: tendency_parse_nml_cpl
  PRIVATE :: compute_eps_and_clear

  ! -------------

CONTAINS

  ! ---------------------------------------------------------------
  ! ===============================================================
  ! MAIN ENTRY POINTS (CALLED FROM MAIN_CONTROL)
  ! ===============================================================
  ! ---------------------------------------------------------------

  SUBROUTINE main_tendency_set_domain(dom_id)

    ! MESSy/BNIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_tracer_mem_bi,    ONLY: ntrac_gp
#ifdef COSMO
    USE data_runcontrol,             ONLY: nnew
    USE data_fields,                 ONLY: pp, t, u, v ,w
#endif
    ! MESSy/SMCL
    USE messy_main_channel_repr,     ONLY: get_representation_id
    USE messy_main_channel,          ONLY: get_channel_object      &
                                         , get_channel_object_info
    USE messy_main_control,          ONLY: main_control_get_context_name

    IMPLICIT NONE

    INTEGER, INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_tendency_set_domain'
    INTEGER                     :: jv, jt, jtstart, jtend
    INTEGER                     :: reprid_mid3d
    CHARACTER(LEN=200)          :: debugstr

    debugstr=substr//main_control_get_context_name()

    !
#ifndef ICON
    CALL get_representation_id(status, 'GP_3D_MID', reprid_mid3d)
#else
    CALL get_representation_id(status, 'UCELL_3D_MID', reprid_mid3d)
#endif
    CALL channel_halt(substr, status)
    !
#ifndef COSMO
    DO jv = NProgVMax, -1
       IF (progvar(jv)%lexists) THEN
          CALL get_channel_object(status                &
             , progvar(jv)%var%cha, progvar(jv)%var%obj &
             , p3=progvar(jv)%varm1ptr, dom_id=dom_id )
          CALL channel_halt(debugstr, status)
          CALL get_channel_object(status                    &
             , progvar(jv)%varte%cha, progvar(jv)%varte%obj &
             , p3=progvar(jv)%varteptr, dom_id=dom_id )
          CALL channel_halt(debugstr, status)
          CALL get_channel_object_info(status                  &
               , progvar(jv)%varte%cha, progvar(jv)%varte%obj &
               , reprid=progvar(jv)%reprid)
          CALL channel_halt(debugstr, status)
       END IF
    END DO

#else
    ! COSMO applies most of the tendencies directly, therefore
    ! we need here a pointer to the actual nnew time-level
    DO jv = NProgVMax, -1
       IF (progvar(jv)%lexists) THEN
          CALL get_channel_object(status                    &
             , progvar(jv)%varte%cha, progvar(jv)%varte%obj &
             , p3=progvar(jv)%varteptr)
          CALL channel_halt(debugstr, status)
          SELECT CASE(jv)
          CASE(mtend_id_pp)
             progvar(jv)%varm1ptr => pp(:,:,:,nnew)
          CASE(mtend_id_t)
             progvar(jv)%varm1ptr => t(:,:,:,nnew)
          CASE(mtend_id_u)
             progvar(jv)%varm1ptr => u(:,:,:,nnew)
          CASE(mtend_id_v)
             progvar(jv)%varm1ptr => v(:,:,:,nnew)
          CASE(mtend_id_w)
             progvar(jv)%varm1ptr => w(:,:,:,nnew)
          END SELECT
       END IF
    END DO

#endif

    ! SET PROGVAR POINTERS
    IF (lany_tracer(dom_curid)) THEN
       DO jt = 1, ntrac_gp
          progvar(jt)%lexists  = .TRUE.
          progvar(jt)%varm1ptr => xtm1(_RI_XYZN_(:,:,:,jt))
          progvar(jt)%varteptr => xtte(_RI_XYZN_(:,:,:,jt))
          progvar(jt)%reprid   = reprid_mid3d
       END DO
       jtstart = ntrac_gp+1
    ELSE
       jtstart = 1
    END IF
    jtend = UBOUND(progvar,1)
    DO jt = jtstart, jtend
       progvar(jt)%lexists  = .FALSE.
       NULLIFY(progvar(jt)%varm1ptr)
       NULLIFY(progvar(jt)%varteptr)
    END DO

  END SUBROUTINE main_tendency_set_domain
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------
  SUBROUTINE main_tendency_initialize(flag)

    USE messy_main_mpi_bi,         ONLY: p_parallel_io, p_bcast
    USE messy_main_tools,          ONLY: find_next_free_unit   &
                                       , domains_from_string
#ifdef COSMO
    USE data_modelconfig,          ONLY: idt_qv, idt_qi, idt_qc &
                                       , idt_qr, idt_qs, idt_qg
    USE data_runcontrol,           ONLY: itype_gscp, ltur, lgsp  &
                                       , itype_conv, itype_turb  &
                                       , l3dturb, lhordiff, lsso, itype_vdif
#ifdef COSMOv5s5
    USE data_runcontrol,           ONLY: loutput_diab
#endif
#endif

    IMPLICIT NONE
    INTEGER, INTENT(IN)         :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='tendency_initialize'
    INTEGER                     :: status
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: jd
    INTEGER                     :: tmp_handle

    INTEGER                     :: num, id, jh
    CHARACTER(LEN=2)            :: hn='00'
    CHARACTER(LEN=32)           :: handlename=''
    CHARACTER(LEN=32)           :: varname
    CHARACTER(LEN=3),  DIMENSION(:), POINTER :: domname => NULL()
    INTEGER,           DIMENSION(:), POINTER :: domnum  => NULL()

    CALL start_message_bi(modstr,'INITIALIZATION',substr)

    SELECT CASE(flag)

    CASE(1)
       ! INITIALIZE CPL
       IF (p_parallel_io) THEN
          iou = find_next_free_unit(100,200)
          ! *** CALL TENDENCY CORE ROUTINE:
          CALL tendency_read_nml_cpl(status, iou)
          IF (status /= 0) CALL error_bi(' ', substr)
       END IF

       DO id = 1, n_dom
          CALL p_bcast(l_closure(id),   p_io)
          CALL p_bcast(l_clos_diag(id), p_io)
          CALL p_bcast(l_full_diag(id), p_io)
       END DO
    !-----------

    ! COUNT VALID ENTRIES ON I/O PE ...
    IF (p_parallel_io) THEN
       N_TDIAG = 1
       DO jd=1,NMAX_TDIAG
          IF (TDIAG(jd)%variable_name == '') CYCLE
          CALL domains_from_string(status, TDIAG(jd)%variable_name,n_dom,num)
          IF (status /= 0) CALL error_bi('error in TDIAG namelist', substr)
          N_TDIAG = N_TDIAG + num
       END DO
       N_TDIAG = N_TDIAG - 1 !=> number of variables in nml
    ENDIF
    ! ... AND BROADCAST RESULT
    CALL p_bcast(N_TDIAG, p_io)
    ! ALLOCATE WORK SPACE
    ALLOCATE(XTDIAG(N_TDIAG))

    ! COPY ENTRIES TO WORK SPACE VARIABLE
    IF (p_parallel_io) THEN
       N_TDIAG = 1
       DO jd=1,NMAX_TDIAG
          IF (TDIAG(jd)%variable_name == '') CYCLE

          CALL domains_from_string(status, TDIAG(jd)%variable_name,n_dom,num &
               ,varname,domname,domnum)
          IF (status /= 0) CALL error_bi('error in TDIAG namelist', substr)
          DO id = 1,SIZE(domname)
             XTDIAG(N_TDIAG)%io%variable_name = TRIM(ADJUSTL(varname))
             XTDIAG(N_TDIAG)%io%domain_name   = TRIM(domname(id))
             XTDIAG(N_TDIAG)%io%domain_idx    = domnum(id)
             XTDIAG(N_TDIAG)%io%process_name = &
                  TRIM(ADJUSTL(TDIAG(jd)%process_name))
             IF (n_dom == 1) THEN
                WRITE(*,*) N_TDIAG, TRIM(XTDIAG(N_TDIAG)%io%variable_name), ':'&
                     ,TRIM(XTDIAG(N_TDIAG)%io%process_name)
             ELSE
                WRITE(*,*) N_TDIAG, TRIM(XTDIAG(N_TDIAG)%io%variable_name), '-'&
                     , TRIM(XTDIAG(N_TDIAG)%io%domain_name),':' &
                     ,TRIM(XTDIAG(N_TDIAG)%io%process_name)
             END IF
             N_TDIAG = N_TDIAG + 1
          END DO

          DEALLOCATE(domname, domnum)
          NULLIFY(domname)
          NULLIFY(domnum)
       END DO
       N_TDIAG = N_TDIAG - 1 !=> number of variables in nml

    ENDIF

    DO jd = 1, N_TDIAG
       CALL p_bcast(XTDIAG(jd)%io%variable_name, p_io)
       CALL p_bcast(XTDIAG(jd)%io%process_name,  p_io)
       CALL p_bcast(XTDIAG(jd)%io%domain_name,   p_io)
       CALL p_bcast(XTDIAG(jd)%io%domain_idx,    p_io)
    ENDDO

 CASE (2)
    ! overwrite humidity mtend_ids for basemodels where those are tracers:
#if defined(COSMO)
    mtend_id_q  = idt_qv
    mtend_id_xl = idt_qc
    mtend_id_xi = idt_qi
#endif

    DO jh=1, NMAX_HANDLE
       ALLOCATE(proc(jh)%p(n_dom))
    END DO

    ! RESERVE 2 "process records" for
    ! - internal summation
    ! - calculation of difference to total internal tendency
    IF(ANY(l_closure)) THEN
       I_HANDLE_SUM  = mtend_get_handle('mtend_sum')
       I_HANDLE_DIFF = mtend_get_handle('mtend_diff')
       DO id= 1, n_dom
          IF (l_closure(id)) THEN
             ! --------- DIFF ------------
             CALL mtend_register(I_HANDLE_DIFF, mtend_id_t,  dom_id=id)
#ifndef ICON
             CALL mtend_register(I_HANDLE_DIFF, mtend_id_u,  dom_id=id)
             CALL mtend_register(I_HANDLE_DIFF, mtend_id_v,  dom_id=id)
#else
             ! is it possible for ICON to directly map mtend_id_q to tracer index as for COSMO ?
             CALL mtend_register(I_HANDLE_DIFF, mtend_id_q,  dom_id=id)
#endif
#ifdef COSMO
             CALL mtend_register(I_HANDLE_DIFF, mtend_id_w,  dom_id=id)
             CALL mtend_register(I_HANDLE_DIFF, mtend_id_pp, dom_id=id)
#endif
#ifdef ECHAM5
             CALL mtend_register(I_HANDLE_DIFF,mtend_id_q,   dom_id=id)
             CALL mtend_register(I_HANDLE_DIFF,mtend_id_xl,  dom_id=id)
             CALL mtend_register(I_HANDLE_DIFF,mtend_id_xi,  dom_id=id)
#endif
             CALL mtend_register(I_HANDLE_DIFF, mtend_id_tracer,  dom_id=id)

             ! --------- SUM ------------

             CALL mtend_register(I_HANDLE_SUM, mtend_id_t,  dom_id=id)
#ifndef ICON
             CALL mtend_register(I_HANDLE_SUM, mtend_id_u,  dom_id=id)
             CALL mtend_register(I_HANDLE_SUM, mtend_id_v,  dom_id=id)
#else
             ! is it possible for ICON to directly map mtend_id_q to tracer index as for COSMO ?
             CALL mtend_register(I_HANDLE_SUM, mtend_id_q,  dom_id=id)
#endif
#ifdef COSMO
             CALL mtend_register(I_HANDLE_SUM, mtend_id_w,  dom_id=id)
             CALL mtend_register(I_HANDLE_SUM, mtend_id_pp, dom_id=id)
#endif
#ifdef ECHAM5
             CALL mtend_register(I_HANDLE_SUM,mtend_id_q,   dom_id=id)
             CALL mtend_register(I_HANDLE_SUM,mtend_id_xl,  dom_id=id)
             CALL mtend_register(I_HANDLE_SUM,mtend_id_xi,  dom_id=id)
#endif
             CALL mtend_register(I_HANDLE_SUM, mtend_id_tracer,  dom_id=id)
          END IF
       END DO
    END IF

#ifdef ECHAM5
    CALL start_message_bi(modstr,'REQUESTS FOR ECHAM5',substr)

    ! ADVECTION
    ! NOTE: The call sequence is:
    !   CALL init_coupling (-> call messy_init_coupling)
    !   CALL init_tpcore
    !   Thus, it is to late to set the requests in init_tpcore, which is
    !   the first entry point of the FFSL advection, because the requests
    !   determine the channel objects below.
    ! Solution: set the requests already here!

    !MO_TPCORE AKA ADVECTION
    tmp_handle=mtend_get_handle('advect')  ! see modstr in mo_tpcore
    CALL mtend_register(tmp_handle,mtend_id_q)
    CALL mtend_register(tmp_handle,mtend_id_xl)
    CALL mtend_register(tmp_handle,mtend_id_xi)
    CALL mtend_register(tmp_handle,mtend_id_tracer)

    !DYN
    tmp_handle=mtend_get_handle('dyn') ! see modstr in dyn
    CALL mtend_register (tmp_handle,mtend_id_t)
    CALL mtend_register (tmp_handle,mtend_id_u)
    CALL mtend_register (tmp_handle,mtend_id_v)
    CALL end_message_bi(modstr,'REQUESTS FOR BASEMODEL ECHAM5',substr)
!ECHAM5
#endif

#ifdef COSMO
    CALL start_message_bi(modstr,'REQUESTS FOR BASEMODEL COSMO',substr)
    ! ADD BASEMODEL PROCESSES here
    ! -- total tendency
    I_HANDLE_BMTT=mtend_get_handle('total')
    CALL mtend_register (I_HANDLE_BMTT,mtend_id_t)
    CALL mtend_register (I_HANDLE_BMTT,mtend_id_u)
    CALL mtend_register (I_HANDLE_BMTT,mtend_id_v)
    CALL mtend_register (I_HANDLE_BMTT,mtend_id_w)
    CALL mtend_register (I_HANDLE_BMTT,mtend_id_pp)
    CALL mtend_register (I_HANDLE_BMTT,mtend_id_tracer)

    ! -- src_relaxation
    tmp_handle=mtend_get_handle('damping')
    CALL mtend_register (tmp_handle,mtend_id_t)
    CALL mtend_register (tmp_handle,mtend_id_u)
    CALL mtend_register (tmp_handle,mtend_id_v)
    CALL mtend_register (tmp_handle,mtend_id_w)
    CALL mtend_register (tmp_handle,mtend_id_pp)
    CALL mtend_register (tmp_handle,mtend_id_tracer)

    tmp_handle=mtend_get_handle('lbc')
    CALL mtend_register (tmp_handle,mtend_id_t)
    CALL mtend_register (tmp_handle,mtend_id_u)
    CALL mtend_register (tmp_handle,mtend_id_v)
    CALL mtend_register (tmp_handle,mtend_id_w)
    CALL mtend_register (tmp_handle,mtend_id_pp)
    CALL mtend_register (tmp_handle,mtend_id_tracer)
    tmp_handle=mtend_get_handle('lbczero')
    CALL mtend_register (tmp_handle,mtend_id_t)
    CALL mtend_register (tmp_handle,mtend_id_u)
    CALL mtend_register (tmp_handle,mtend_id_v)
    CALL mtend_register (tmp_handle,mtend_id_w)
    CALL mtend_register (tmp_handle,mtend_id_pp)
    CALL mtend_register (tmp_handle,mtend_id_tracer)
    tmp_handle=mtend_get_handle('nullify_tracers')
    CALL mtend_register (tmp_handle,mtend_id_tracer)

    ! -- advection
    tmp_handle=mtend_get_handle('advect')
    CALL mtend_register (tmp_handle,mtend_id_tracer)
    CALL mtend_register (tmp_handle,mtend_id_u)
    CALL mtend_register (tmp_handle,mtend_id_v)
    CALL mtend_register (tmp_handle,mtend_id_w)
    CALL mtend_register (tmp_handle,mtend_id_t)
    CALL mtend_register (tmp_handle,mtend_id_pp)

    ! -- dynamics (runge_kutta /slow_tendencies)
    tmp_handle=mtend_get_handle('clip_rk')
    CALL mtend_register (tmp_handle,idt_qv)
    CALL mtend_register (tmp_handle,idt_qc)
    CALL mtend_register (tmp_handle,idt_qr)
    IF (itype_gscp >= 2) &
    CALL mtend_register (tmp_handle,idt_qs)
    IF (itype_gscp >= 3) &
    CALL mtend_register (tmp_handle,idt_qi)
    IF (itype_gscp >= 4) &
    CALL mtend_register (tmp_handle,idt_qg)
    ! stochastic physics tendency
    tmp_handle=mtend_get_handle('RK_sppt')
    CALL mtend_register (tmp_handle,mtend_id_t)
    CALL mtend_register (tmp_handle,mtend_id_u)
    CALL mtend_register (tmp_handle,mtend_id_v)
    ! -- horizontal diffusion
    IF ( l3dturb ) THEN
       tmp_handle=mtend_get_handle('explhor_diff')
       CALL mtend_register (tmp_handle,mtend_id_t)
       CALL mtend_register (tmp_handle,mtend_id_u)
       CALL mtend_register (tmp_handle,mtend_id_v)
       CALL mtend_register (tmp_handle,mtend_id_w)
       CALL mtend_register (tmp_handle,mtend_id_tracer)
    END IF
    IF (lhordiff) THEN
       tmp_handle=mtend_get_handle('horizontal_diff')
       CALL mtend_register (tmp_handle,mtend_id_t)
       CALL mtend_register (tmp_handle,mtend_id_u)
       CALL mtend_register (tmp_handle,mtend_id_v)
       CALL mtend_register (tmp_handle,mtend_id_w)
       CALL mtend_register (tmp_handle,mtend_id_tracer)
       CALL mtend_register (tmp_handle,mtend_id_pp)
    END IF
    ! -- vertical diffusion
    IF (itype_vdif == -1) THEN
       tmp_handle=mtend_get_handle('vdiff')
       CALL mtend_register (tmp_handle,mtend_id_tracer)
       tmp_handle=mtend_get_handle('sppt_vdiff')
       CALL mtend_register (tmp_handle,mtend_id_tracer)
    END IF
    tmp_handle=mtend_get_handle('clipslowtend')
    CALL mtend_register (tmp_handle,idt_qv)
    CALL mtend_register (tmp_handle,idt_qc)
    CALL mtend_register (tmp_handle,idt_qr)
    IF (itype_gscp >= 2) &
         CALL mtend_register (tmp_handle,idt_qs)
    IF (itype_gscp >= 3) &
         CALL mtend_register (tmp_handle,idt_qi)
    IF (itype_gscp >= 4) &
         CALL mtend_register (tmp_handle,idt_qg)

    tmp_handle=mtend_get_handle('impl_vdiff')
    CALL mtend_register (tmp_handle,mtend_id_t)
    CALL mtend_register (tmp_handle,mtend_id_u)
    CALL mtend_register (tmp_handle,mtend_id_v)
    CALL mtend_register (tmp_handle,mtend_id_w)

    tmp_handle=mtend_get_handle('vdiff_coriolis')
    CALL mtend_register (tmp_handle,mtend_id_w)

    ! moisture divergence
    !tmp_handle=mtend_get_handle('dqvdt_vdiff')

    !tmp_handle=mtend_get_handle('dqvdt_advect')

   ! tmp_handle=mtend_get_handle('dqvdt_conv')

    ! latent heat  (src_runge_kutta)
    tmp_handle=mtend_get_handle('lh')
    CALL mtend_register (tmp_handle,mtend_id_t)

    ! saturation adjustment
    tmp_handle=mtend_get_handle('satad_relax')
    CALL mtend_register (tmp_handle,idt_qv)
    CALL mtend_register (tmp_handle,idt_qc)
    CALL mtend_register (tmp_handle,mtend_id_t)
    tmp_handle=mtend_get_handle('satad_dyncore')
    CALL mtend_register (tmp_handle,idt_qv)
    CALL mtend_register (tmp_handle,idt_qc)
    CALL mtend_register (tmp_handle,mtend_id_t)

    ! -- microphysics
    tmp_handle=mtend_get_handle('gscp')
    CALL mtend_register (tmp_handle,mtend_id_t)
    CALL mtend_register (tmp_handle,idt_qv)
    CALL mtend_register (tmp_handle,idt_qc)
    CALL mtend_register (tmp_handle,idt_qr)
    IF (itype_gscp >= 2) &
    CALL mtend_register (tmp_handle,idt_qs)
    IF (itype_gscp >= 3) &
    CALL mtend_register (tmp_handle,idt_qi)
    IF (itype_gscp >= 4) &
    CALL mtend_register (tmp_handle,idt_qg)

    ! -- convection
    tmp_handle=mtend_get_handle('conv')
    IF (itype_conv==0 .OR. itype_conv == 3) THEN
       CALL mtend_register (tmp_handle,mtend_id_t)
       CALL mtend_register (tmp_handle,idt_qv)
       CALL mtend_register (tmp_handle,mtend_id_u)
       CALL mtend_register (tmp_handle,mtend_id_v)
       !CALL mtend_register (tmp_handle,mtend_id_tke)
       !IF (itype_conv == 0) THEN
       CALL mtend_register (tmp_handle,idt_qc)
       CALL mtend_register (tmp_handle,idt_qi)
       !END IF
#ifdef COSMOv5s5
       IF (loutput_diab) THEN
          tmp_handle=mtend_get_handle('conv_diabatic')
          CALL mtend_register (tmp_handle,mtend_id_t)
       END IF
#endif
    END IF
    ! -- turbulence
    IF (itype_turb==3) THEN
       tmp_handle=mtend_get_handle('turbdiff')
       CALL mtend_register (tmp_handle,mtend_id_u)
       CALL mtend_register (tmp_handle,mtend_id_v)
       CALL mtend_register (tmp_handle,mtend_id_t)
       !CALL mtend_register (tmp_handle,idt_qv)
       !CALL mtend_register (tmp_handle,idt_qc)
    END IF
    IF (itype_turb==3 .AND. itype_vdif==1) THEN
       tmp_handle=mtend_get_handle('vertdiff')
       CALL mtend_register (tmp_handle,mtend_id_u)
       CALL mtend_register (tmp_handle,mtend_id_v)
       CALL mtend_register (tmp_handle,mtend_id_t)
       CALL mtend_register (tmp_handle,mtend_id_tracer)
    END IF
    ! -- TERRA (SOIL)
    ! -- lake
    ! -- SSO
    IF (lsso) THEN
       tmp_handle=mtend_get_handle('sso')
       CALL mtend_register (tmp_handle,mtend_id_t)
       CALL mtend_register (tmp_handle,mtend_id_u)
       CALL mtend_register (tmp_handle,mtend_id_v)
    ENDIF
    ! -- radiation
    tmp_handle=mtend_get_handle('radiation')
    CALL mtend_register (tmp_handle,mtend_id_t)
    ! --  ...
    CALL end_message_bi(modstr,'REQUESTS FOR BASEMODEL COSMO',substr)
#endif

#ifdef ICON
    CALL start_message_bi(modstr,'REQUESTS FOR BASEMODEL ICON',substr)
    ! -- convection

    CALL end_message_bi(modstr,'REQUESTS FOR BASEMODEL ICON',substr)
#endif

 END SELECT

    CALL end_message_bi(modstr,'INITIALIZATION',substr)

    RETURN
  END SUBROUTINE main_tendency_initialize
  ! ---------------------------------------------------------------

#ifdef COSMO
  ! ---------------------------------------------------------------
  SUBROUTINE main_tendency_new_tracer

    IMPLICIT NONE

    CALL main_tendency_initialize(2)

  END SUBROUTINE main_tendency_new_tracer
  ! ---------------------------------------------------------------
#endif

  ! ---------------------------------------------------------------
  SUBROUTINE main_tendency_init_coupling

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_tracer_mem_bi,    ONLY: ntrac_gp, GPTRSTR
#ifdef ICON
    USE messy_main_tracer_mem_bi,    ONLY: l_ntrac_icon, L_GPTRSTR
    USE messy_main_tracer_tools_bi,  ONLY: main_tracer_set_domain
    USE messy_main_channel_bi,       ONLY: main_channel_set_domain
#endif
    USE messy_main_constants_mem,    ONLY: STRLEN_MEDIUM
    USE messy_main_tools,            ONLY: int2str
    USE messy_main_tracer,           ONLY: get_tracer, STRLEN_TRSET
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt

    USE messy_main_channel_repr,     ONLY: get_representation_id
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_channel_object_reference    &
                                         , new_attribute
    USE messy_main_mpi_bi,           ONLY: p_parallel_io


    implicit none

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER      :: substr = 'main_tendency_init_coupling'
    INTEGER                          :: jh, jv, jt, status
    CHARACTER(LEN=2*STRLEN_MEDIUM+1) :: name    = '' ! name_subname (tracer)
    CHARACTER(LEN=STRLEN_MEDIUM)     :: unit     = ''
    CHARACTER(LEN=STRLEN_MEDIUM)     :: unit1    = ''
    INTEGER                          :: jd, jc, ih, ic, id, ix
    LOGICAL                          :: lunacc
    LOGICAL                          :: lvarreq,ltrareq
    LOGICAL, DIMENSION(NProgVmax:-1) :: lvmod_loc = .FALSE.
    CHARACTER(LEN=STRLEN_TRSET)      :: GP_TRSTR
    INTEGER                          :: ntrac_max = 0
    CHARACTER(LEN=2*STRLEN_MEDIUM+1) :: objname    = '' ! name_subname (tracer)
    CHARACTER(len=3)                 :: int_string
    INTEGER                          :: reprid_mid3d

    CALL start_message_bi(modstr,'INIT COUPLING',substr)
    DO id = 1, n_dom
       ! get representation id
#ifndef ICON
       CALL get_representation_id(status, 'GP_3D_MID', reprid_mid3d)
#else
       CALL get_representation_id(status, 'UCELL_3D_MID' &
            , reprid_mid3d, dom_id=id)
#endif
       CALL channel_halt(substr, status)

       CALL new_channel(status, modstr, reprid=reprid_mid3d &
#ifdef ICON
            , dom_id=id &
#endif
            )
       CALL channel_halt(substr//'1', status)

       DO ic = 1, NREQ
          CALL new_channel(status, modstr//chpf(ic), reprid=reprid_mid3d &
#ifdef ICON
            , dom_id=id &
#endif
            )
          CALL channel_halt(substr//'2', status)
       END DO
    ENDDO

    !add loops
    ntrac_max=0
    ALLOCATE(lany_tracer(n_dom))
    lany_tracer = .FALSE.
#ifndef ICON
    DO jh = 1, N_HANDLE
       IF (proc(jh)%p(1)%any_tracer) THEN
          proc(jh)%p(1)%ntrac_proc=ntrac_gp
          ntrac_max = ntrac_gp
          lany_tracer(:) = .TRUE.
       END IF
    END DO
#else
    DO id = 1, n_dom
       CALL main_tracer_set_domain(id)
       DO jh = 1, N_HANDLE
          IF (proc(jh)%p(id)%any_tracer) THEN
             proc(jh)%p(id)%ntrac_proc=ntrac_gp
             ntrac_max=MAX(ntrac_max, ntrac_gp)
             lany_tracer(id) = .TRUE.
          END IF
       END DO
    END DO
#endif

    ALLOCATE(progvar(NProgVMax:ntrac_max))
    ! SET POINTER FOR TRACERs, if any tracers requested
    ! moved to set_patch called in set_basemodels_progvarset for the first time

    ! SET BASEMODEL SPECIFIC PROGNOSTIC VARIABLES
    CALL  set_basemodels_progvarset

    !call to read tendency nml and install channels for sums
    CALL tendency_parse_nml_cpl

    handle_loop: DO jh = 1, N_HANDLE
       IF (p_parallel_io) &
            WRITE(*,*) 'PROCESS ',TRIM(proc(jh)%name), ' (HANDLE=',jh,') ...'

       patch_loop: DO id = 1, n_dom
#ifdef ICON
          CALL main_channel_set_domain(id)
          CALL main_tracer_set_domain(id)
          CALL main_tendency_set_domain(id)
#endif
          ! ALLOCATE POINTER ARRAY AND NULLIFY 3D POINTERS FOR VARIABLES
          DO ic = 1, NCASES
             ALLOCATE(proc(jh)%p(id)%case(ic)%pvar(NProgVmax:proc(jh)%p(id)%ntrac_proc))
             DO jv= NProgVmax, proc(jh)%p(id)%ntrac_proc
                NULLIFY(proc(jh)%p(id)%case(ic)%pvar(jv)%ptr)
             END DO
             ALLOCATE(proc(jh)%p(id)%case(ic)%lvar(NProgVmax:proc(jh)%p(id)%ntrac_proc))
             proc(jh)%p(id)%case(ic)%lvar(:) = .FALSE.
             ALLOCATE(proc(jh)%p(id)%case(ic)%lextern(NProgVmax:proc(jh)%p(id)%ntrac_proc))
             proc(jh)%p(id)%case(ic)%lextern(:) = .FALSE.
          ENDDO

          ! GET TRACERS, ALLOCATE POINTER ARRAY AND NULLIFY 3D POINTERS FOR
          ! TRACERS
          any_tracer: IF (proc(jh)%p(id)%any_tracer) THEN

             lvmod_loc = proc(jh)%p(id)%lvmod
             DEALLOCATE( proc(jh)%p(id)%lvmod)

             ALLOCATE(proc(jh)%p(id)%lvmod(NProgVmax:proc(jh)%p(id)%ntrac_proc))
             proc(jh)%p(id)%lvmod(:) = .FALSE.
             proc(jh)%p(id)%lvmod(NProgVmax:-1) =  lvmod_loc
             lvmod_loc = .FALSE.

             IF (SIZE(proc(jh)%p(id)%idx) > 0) THEN
                ! only selected tracers
                tracer_loop: DO jt=1, SIZE(proc(jh)%p(id)%idx)
                   proc(jh)%p(id)%lvmod(proc(jh)%p(id)%idx(jt)) = .TRUE.
                END DO tracer_loop
             ELSE
                ! all tracers
                proc(jh)%p(id)%lvmod(1:proc(jh)%p(id)%ntrac_proc) = .TRUE.
             END IF
             !=============
          ENDIF any_tracer

          ! --------------------------------------------------
          ! INTERPRET CPL NAMELIST
          ! --------------------------------------------------
          DO ic = 1, NREQ

             SELECT CASE(ic)

             CASE(IR_NML)
                !-------------------------------------------------------------
                ! case for specific (sums of) tendencies requested by namelist
                !-------------------------------------------------------------

                variable_loop_nml: DO jv = NProgVmax, proc(jh)%p(id)%ntrac_proc

                   !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
                   !loop structure from subroutine parse_nml_cpl
                   !to compare handles and point on channel of sum for variables
                   IF (jv == 0) CYCLE
                   IF (jv > 0) THEN

                      CALL get_tracer(status, GPTRSTR, jv, fullname=name, &
                           unit=unit1)
                      CALL tracer_halt(substr, status)
                   ELSE
                      name = TRIM(pvarname(jv))
                   END IF


                   valid_var_nml: IF (proc(jh)%p(id)%lvmod(jv)) THEN

                      IF (p_parallel_io) &
                        WRITE(*,*) ' ... (IC_NML) modifies variable ',TRIM(name)

                      !parse loops
                      lvarreq = .FALSE.
                      tdiag_loop_var: DO jd=1,N_TDIAG
                         !-------------
                         ! skip internal special accounting handles ...
                         IF(l_closure(id)) THEN
                            IF (jh == I_HANDLE_SUM)  CYCLE
                            IF (jh == I_HANDLE_DIFF) CYCLE
                         ENDIF
                         !-------------
                         !Comparing variable names fit
                         var_name: IF (TRIM(name) .eq.               &
                              TRIM(XTDIAG(jd)%io%variable_name)      &
                              .AND. id == XTDIAG(jd)%io%domain_idx   &
                              ) THEN

                            lvarreq = .TRUE.
                            IF (p_parallel_io) &
                                WRITE(*,*) ' ... ... tendency diagnosed in ... '

                            ! POINT TO UNACCOUNTED FIRST ...
                            proc(jh)%p(id)%case(IC_NML)%pvar(jv)%ptr => &
                                 XTDIAG(jd)%obj(XTDIAG(jd)%ntclass)%ptr(:,:,:)
                            lunacc = .TRUE.

                            ! .. THEN LOOP OVER ALL CLASSES ...
                            class_loop:  DO jc = 1, XTDIAG(jd)%ntclass

                               p_handle_loop: DO  ih = 1, XTDIAG(jd)%nproc(jc)
                                  !(=nr2)
                                  ! ... AND CHECK FOR MATCHING HANDLE WITHIN CLASS

                                  IF (proc(jh)%handle .eq. &
                                       XTDIAG(jd)%iproc(jc)%ptr(ih) .AND. &
                                       id == XTDIAG(jd)%io%domain_idx) THEN
                                     IF (p_parallel_io) THEN
                                        WRITE(*,*) ' ... ... ... user request ',&
                                             jd &
                                             ,'; class ',jc,'; summand ',ih &
                                             ,' (handle=',&
                                             XTDIAG(jd)%iproc(jc)%ptr(ih),')'
                                     END IF
                                     proc(jh)%p(id)%case(IC_NML)%pvar(jv)%ptr => &
                                          XTDIAG(jd)%obj(jc)%ptr(:,:,:)
                                     lunacc = .FALSE.

                                     l1proc: IF ( XTDIAG(jd)%nproc(jc) == 1) THEN
                                        ! USE MEMORY FOR TENDENCY
                                        ! BUT ONLY IF XTDIAG FOR ONLY 1 PROCESS
                                        proc(jh)%p(id)%case(IC_GEN)%lvar(jv) = .TRUE.
                                        proc(jh)%p(id)%case(IC_GEN)%pvar(jv)%ptr => &
                                              XTDIAG(jd)%obj(jc)%ptr(:,:,:)

                                        CALL int2str(int_string,jc)
                                        objname=&
                                           TRIM(XTDIAG(jd)%io%variable_name)//'_nml_sum_nr'//TRIM(int_string)
                                        CALL new_channel_object_reference(status,          &
                                             modstr//chpf(IR_NML), TRIM(objname),          &
                                             modstr, TRIM(name)//'_'//TRIM(proc(jh)%name), &
                                             lcopyatt=.TRUE.                               &

#ifdef ICON
                                             , dom_id1=id , dom_id2=id &
#endif
                                             )
                                        CALL channel_halt(substr//'3', status)
                                        proc(jh)%p(id)%case(IC_NML)%lvar(jv) = .FALSE.
                                     ELSE
                                        !
                                        ! This switch indicates that the IC_NML
                                        ! pointer exists, but that the memory is
                                        ! used for several processes at once
                                        proc(jh)%p(id)%case(IC_NML)%lvar(jv) = .TRUE.
                                     END IF l1proc
                                  ENDIF

                               ENDDO p_handle_loop

                            ENDDO class_loop

                            IF (lunacc) THEN
                               IF (p_parallel_io) THEN
                                  WRITE(*,*) ' ... ... ... user request ',jd,&
                                       ' (unaccounted)'
                               ENDIF
                               proc(jh)%p(id)%case(IC_NML)%lvar(jv) = .TRUE.
                            ENDIF

                         END IF var_name

                      ENDDO tdiag_loop_var

                      IF (.NOT. lvarreq) THEN
                         IF (p_parallel_io) THEN
                            WRITE(*,*) ' ... ... tendency not diagnosed'
                         ENDIF
                      END IF

                   ELSE
                      ! NULL-POINTER FOR invalid variables
                      NULLIFY(proc(jh)%p(id)%case(IC_NML)%pvar(jv)%ptr)
                   END IF valid_var_nml
                   !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

                ENDDO variable_loop_nml

             CASE(IR_FUL)
                !---------------------------------------------------------------
                ! case for all tendencies are needed seperately (l_full_diag)
                !--------------------------------------------------------------
                IF (l_full_diag(id) .AND. jh /= I_HANDLE_SUM .AND. &
                     jh /= I_HANDLE_DIFF) THEN

                   variable_loop_ful: DO jv = NProgVmax, proc(jh)%p(id)%ntrac_proc
                      IF (jv == 0) CYCLE
                      IF (jv > 0) THEN
                         CALL get_tracer(status, GPTRSTR, jv, fullname=name, &
                              unit=unit1)
                         CALL tracer_halt(substr, status)
                         unit = '('//TRIM(unit1)//') / s'
                      ELSE
                         name = pvarname(jv)
                         unit = pvarunit(jv)
                      END IF
                      valid_var_ful: IF (proc(jh)%p(id)%lvmod(jv)) THEN
                         ! ORIGINAL: create channel objects for all valid process
                         !           / variable pairs
                         !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                         IF (p_parallel_io) THEN
                         IF (p_parallel_io) &
                              WRITE(*,*) ' ... (IR_FUL) modifies variable ',TRIM(name)
                         END IF
                         lasso: IF (.NOT. proc(jh)%p(id)%case(IC_GEN)%lvar(jv)) THEN
                            CALL new_channel_object(status, modstr          &
                                 , TRIM(name)//'_'//TRIM(proc(jh)%name)     &
                                 , p3 = proc(jh)%p(id)%case(IC_GEN)%pvar(jv)%ptr    &
                                 , reprid=progvar(jv)%reprid                &
                                 )
                            proc(jh)%p(id)%case(IC_GEN)%lvar(jv) = .TRUE.
                            CALL channel_halt(substr//'4', status)
                            CALL new_attribute(status,modstr              ,&
                                 TRIM(name)//'_'//TRIM(proc(jh)%name),     &
                                 'long_name',                              &
                                 c=TRIM(name)//': tendencies: '//TRIM(proc(jh)%name)  &
                                 )
                            CALL channel_halt(substr//'5', status)
                            CALL new_attribute(status,modstr,              &
                                 TRIM(name)//'_'//TRIM(proc(jh)%name),     &
                                 'units',c=TRIM(unit)                      &
                                 )
                            CALL channel_halt(substr//'6', status)
                         END IF lasso

                         CALL new_channel_object_reference(status,          &
                              modstr, TRIM(name)//'_'//TRIM(proc(jh)%name), &
                              modstr//chpf(IR_FUL),                         &
                              TRIM(name)//'_'//TRIM(proc(jh)%name),         &
                              lcopyatt=.TRUE.                               &
                              )
                         CALL channel_halt(substr//'7', status)

                      ENDIF valid_var_ful
                   ENDDO variable_loop_ful

                   !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
                ENDIF

             CASE(IR_EXC)
                !---------------------------------------------------------------
                ! case for tendencies requested by other submodels via call
                ! mtend_request
                !---------------------------------------------------------------

                ! nothing to do here
                ! everything is done in the subroutine mtend_request

             CASE(IR_CLS)
                !-----------------------------------------------------------------
                ! case for the two special objects i_handle_sum/-diff (l_closure)
                !-----------------------------------------------------------------
                IF(l_closure(id) .AND. (jh == I_HANDLE_SUM .OR. &
                     jh == I_HANDLE_DIFF .OR. jh == I_HANDLE_BMTT)) THEN

                   ! create channel objects for i_handle_sum /-diff with all
                   ! possible variables
                   variable_loop_cls: DO jv = NProgVmax, proc(jh)%p(id)%ntrac_proc

                      IF (jv == 0) CYCLE
                      IF (jv > 0) THEN
                         CALL get_tracer(status, GPTRSTR, jv, fullname=name, &
                              unit=unit1)
                         CALL tracer_halt(substr, status)
                         unit = '('//TRIM(unit1)//') / s'
                      ELSE
                         name = pvarname(jv)
                         unit = pvarunit(jv)
                      END IF
                      valid_var_cls: IF (proc(jh)%p(id)%lvmod(jv)) THEN
                         IF (.NOT. proc(jh)%p(id)%case(IC_GEN)%lvar(jv)) THEN

                            CALL new_channel_object(status, modstr        &
                                 , TRIM(name)//'_'//TRIM(proc(jh)%name)   &
                                 , p3=proc(jh)%p(id)%case(IC_GEN)%pvar(jv)%ptr &
                                 , reprid=progvar(jv)%reprid              &
                                 )
                            CALL channel_halt(substr//'8', status)
                            CALL new_attribute(status,modstr,             &
                                 TRIM(name)//'_'//TRIM(proc(jh)%name),    &
                                 'long_name', c=TRIM(name)//'_te'         &
                                 )
                            CALL channel_halt(substr//'9', status)
                            CALL new_attribute(status,modstr,              &
                                 TRIM(name)//'_'//TRIM(proc(jh)%name),     &
                                 'units',c=TRIM(unit)                      &
                                 )
                            CALL channel_halt(substr//'10', status)
                            proc(jh)%p(id)%case(IC_GEN)%lvar(jv) = .TRUE.
                         END IF
                         CALL new_channel_object_reference(status,          &
                              modstr, TRIM(name)//'_'//TRIM(proc(jh)%name), &
                              modstr//chpf(IR_CLS),                         &
                              TRIM(name)//'_'//TRIM(proc(jh)%name),         &
                              lcopyatt=.TRUE.                               &
                              )
                         CALL channel_halt(substr//'11', status)
                      ENDIF valid_var_cls
                   ENDDO variable_loop_cls

                   !--------------------------

                ENDIF

             END SELECT
          ENDDO
       END DO patch_loop
    END DO handle_loop

    CALL end_message_bi(modstr,'INIT COUPLING',substr)

  END SUBROUTINE main_tendency_init_coupling
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------
  SUBROUTINE main_tendency_global_end

    USE messy_main_grid_def_mem_bi, ONLY: nproma, ngpblks
    USE messy_main_channel_mem,     ONLY: dom_curid
    USE messy_main_timer,           ONLY: lfirst_cycle
    USE messy_main_tracer_mem_bi,   ONLY: ti_gp

    IMPLICIT NONE
    INTEGER            :: jd, jc, jt, jh, jk, n, m, ic, jv
    REAL(kind=dp), DIMENSION(:,:), POINTER :: zzz => NULL()

#ifdef ICON
    ! add all MESSy tracer tendencies to actual tracers
    !TODO CALL icon_tracer_integration
#endif

     IF(lfirst_cycle) THEN
       IF(l_closure(dom_curid)) THEN

          write(6,*) ' '
          write(6,*) '################ Tendency Budget #######################################'
          write(6,*) ' '
          write(6,*) ' ==========> Tendency records have been defined for the following variables and processes: '
          write(6,*) ' '

          DO jd = 1,N_TDIAG
             IF (XTDIAG(jd)%io%domain_idx /= dom_curid) CYCLE
             write(6,*) '  Variables: ', TRIM(XTDIAG(jd)%io%variable_name)
             IF (TRIM(XTDIAG(jd)%io%process_name) == '') THEN
                write(6,*) '  Processes: ', 'All in unaccounted'
             ELSE
                write(6,*) '  Processes: ', TRIM(XTDIAG(jd)%io%process_name)
             ENDIF
          ENDDO

          write(6,*) ' '
          write(6,*) '########################################################################'
          write(6,*) ' '

       ENDIF
    ENDIF

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Calculate difference between total and internal tendencies
    IF(l_closure(dom_curid)) THEN
       !Here: t,q,xl,xi. (u and v are below due to scaling)
       DO jv = NProgVMax, -3
          IF (progvar(jv)%lexists) THEN
             proc(I_HANDLE_DIFF)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(:,:,:)  =  &
                  proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(:,:,:) &
!#if defined(COSMO) || defined(ICON)
#if defined(COSMO)
                  - proc(I_HANDLE_BMTT)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(:,:,:)
#else
                  - progvar(jv)%varteptr(:,:,:)
#endif
             call compute_eps_and_clear (TRIM(progvar(jv)%var%obj)  &
!#if defined(COSMO) || defined(ICON)
#if defined(COSMO)

                  ,  proc(I_HANDLE_BMTT)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(:,:,:), &
#else
                  ,  progvar(jv)%varteptr(:,:,:) ,  &
#endif
                  proc(I_HANDLE_DIFF)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(:,:,:))
          END IF
       END DO

       ! U,v
       ALLOCATE(zzz(nproma,ngpblks))
       IF(wind_is_scaled) THEN
#ifdef WINDSCALING
          zzz = sqcst_2d
#else
          zzz = 1._dp
#endif
          where(zzz(:,:) < 1.0E-40_dp)
             zzz(:,:) = 1._dp
          end where
       ELSE
          zzz(:,:) = 1._dp
       END IF

       DO jv = -2,-1
          IF (progvar(jv)%lexists) THEN
             DO jk=1,nlev
                proc(I_HANDLE_DIFF)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(_RI_XYZ__(:,:,jk)) =&
                     proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(_RI_XYZ__(:,:,jk)) &
!#if !defined(COSMO) && !defined(ICON)
#if !defined(COSMO)
                     - progvar(jv)%varteptr(_RI_XYZ__(:,:,jk))/zzz(:,:)
#else
                - proc(I_HANDLE_BMTT)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr( _RI_XYZ__(:,:,jk) )
#endif
             END DO
             call compute_eps_and_clear (TRIM(progvar(jv)%var%obj)  &
!#if !defined(COSMO) && !defined(ICON)
#if !defined(COSMO)
                  ,  progvar(jv)%varteptr(:,:,:) ,  &
#else
                  , proc(I_HANDLE_BMTT)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(:,:,:), &
#endif
                  proc(I_HANDLE_DIFF)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(:,:,:))
          END IF
       END DO

       ! TRACER, if requested
       DO jv =1, proc(I_HANDLE_DIFF)%p(dom_curid)%ntrac_proc
          IF (progvar(jv)%lexists) THEN
!#if  !defined(COSMO) && !defined(ICON)
#if  !defined(COSMO)
             proc(I_HANDLE_DIFF)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(:,:,:)  =  &
                  proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(:,:,:) &
                  - progvar(jv)%varteptr(:,:,:)
             call compute_eps_and_clear (&
                  '_xtte: '//TRIM(ti_gp(jv)%tp%ident%fullname)  &
                  ,  progvar(jv)%varteptr(:,:,:) ,            &
                  proc(I_HANDLE_DIFF)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(:,:,:))
#else
             proc(I_HANDLE_DIFF)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(:,:,:)  =  &
                  proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(:,:,:) &
                  - proc(I_HANDLE_BMTT)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(:,:,:)
             call compute_eps_and_clear (&
                  '_xtte: '//TRIM(ti_gp(jv)%tp%ident%fullname)                  &
                  , proc(I_HANDLE_BMTT)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(:,:,:) &
                  , proc(I_HANDLE_DIFF)%p(dom_curid)%case(IC_GEN)%pvar(jv)%ptr(:,:,:) )
#endif

          END IF
       END DO

       DEALLOCATE(zzz)
       NULLIFY(zzz)
    ENDIF

    RETURN

  END SUBROUTINE main_tendency_global_end
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------
  SUBROUTINE main_tendency_reset(dom_id)

    USE messy_main_channel_mem, ONLY: dom_curid

    IMPLICIT NONE

    INTEGER, INTENT(in), OPTIONAL :: dom_id

    integer                    :: jh, jv, jt, ic, pix
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tendency_reset'

    IF (PRESENT(dom_id)) THEN
       pix = dom_id
    ELSE
       pix = dom_curid
    END IF

    process_loop: DO jh = 1, N_HANDLE

       variable_loop: DO jv=NProgVmax,proc(jh)%p(pix)%ntrac_proc
          IF(proc(jh)%p(pix)%lvmod(jv)) THEN
             ! non-tracer prognostic variables
             DO ic = 1, NCASES
                IF (proc(jh)%p(pix)%case(ic)%lvar(jv)) &
                     proc(jh)%p(pix)%case(ic)%pvar(jv)%ptr(:,:,:) = 0.0_dp
             END DO
          ENDIF
       ENDDO variable_loop

    END DO process_loop

#ifdef COSMO
    DO jv=NProgVmax,-1
       IF (progvar(jv)%lexists)  progvar(jv)%varteptr = 0._dp
    END DO
#endif
#ifdef ICON
    ! reset MESSy tracer tendencies for current patch
    xtte(:,:,:,:) = 0.0_dp
#endif

  END SUBROUTINE main_tendency_reset
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------
  SUBROUTINE main_tendency_free_memory

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tendency_free_memory'
    INTEGER                     :: jh, jd, jc, ic, id

    handle_loop: DO jh = 1, N_HANDLE

       patch_loop: DO id = 1, n_dom

          any_tracer: IF (proc(jh)%p(id)%any_tracer) THEN

             IF (ASSOCIATED(proc(jh)%p(id)%idx)) THEN
                DEALLOCATE(proc(jh)%p(id)%idx)
                NULLIFY(proc(jh)%p(id)%idx)
             END IF

          END IF any_tracer

          IF (ASSOCIATED(proc(jh)%p(id)%lvmod)) THEN
             DEALLOCATE(proc(jh)%p(id)%lvmod)
             NULLIFY(proc(jh)%p(id)%lvmod)
          END IF

          DO ic = 1, NCASES
             IF (ASSOCIATED(proc(jh)%p(id)%case(ic)%pvar)) THEN
                DEALLOCATE(proc(jh)%p(id)%case(ic)%pvar)
                NULLIFY(proc(jh)%p(id)%case(ic)%pvar)
             END IF

             IF (ASSOCIATED(proc(jh)%p(id)%case(ic)%lvar)) THEN
                DEALLOCATE(proc(jh)%p(id)%case(ic)%lvar)
                NULLIFY(proc(jh)%p(id)%case(ic)%lvar)
             END IF

          ENDDO
      END DO patch_loop
    END DO handle_loop

    DO jd=1, SIZE(XTDIAG)
       IF (ASSOCIATED(XTDIAG(jd)%nproc)) THEN
          DEALLOCATE(XTDIAG(jd)%nproc)
          NULLIFY(XTDIAG(jd)%nproc)
       END IF
       IF (ASSOCIATED(XTDIAG(jd)%obj)) THEN
          DEALLOCATE(XTDIAG(jd)%obj)
          NULLIFY(XTDIAG(jd)%obj)
       END IF
       IF (ASSOCIATED(XTDIAG(jd)%iproc)) THEN
          DO jc=1, SIZE(XTDIAG(jd)%iproc)
             IF (ASSOCIATED(XTDIAG(jd)%iproc(jc)%ptr)) THEN
                DEALLOCATE(XTDIAG(jd)%iproc(jc)%ptr)
                NULLIFY(XTDIAG(jd)%iproc(jc)%ptr)
             ENDIF
          END DO
          DEALLOCATE(XTDIAG(jd)%iproc)
          NULLIFY(XTDIAG(jd)%iproc)
       END IF
    END DO

    IF (ASSOCIATED(XTDIAG)) THEN
       DEALLOCATE(XTDIAG)
       NULLIFY(XTDIAG)
    END IF

    DEALLOCATE(progvar)
    NULLIFY(PROGVAR)

    DEALLOCATE(lany_tracer)

  END SUBROUTINE main_tendency_free_memory
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------
  ! ===============================================================
  ! ENTRY POINTS CALLED FROM SMIL
  ! ===============================================================
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------

  INTEGER FUNCTION mtend_get_handle(name, lnew)

    ! THIS FUNCTION MUST BE CALLED TO "REGISTER" A PROCESS (SUBMODEL)
    ! IN ORDER TO ALLOW TENDENCY CHANGES OF PROGNOSTIC VARIABLES AND
    ! TRACERS

    IMPLICIT NONE
    INTRINSIC :: TRIM, ADJUSTL

    ! I/O
    CHARACTER(LEN=*), INTENT(IN) :: name ! name of the calling submodel
    LOGICAL,          INTENT(IN), OPTIONAL :: lnew ! allow new handle?

    ! which tendencies should be modified
    ! LOCAL
    INTEGER                               :: jh,ind, jv, id
    LOGICAL                               :: NotFound
    LOGICAL                               :: zlnew
    CHARACTER(LEN=*),PARAMETER            :: substr = 'mtend_get_handle'

    IF (PRESENT(lnew)) THEN
       zlnew = lnew
    ELSE
       zlnew = .TRUE.  ! default
    END IF

    NotFound = .true.
    do jh=1,N_HANDLE
       if(TRIM(ADJUSTL(name)) == TRIM(ADJUSTL(proc(jh)%name))) then
          ind      = proc(jh)%handle
          NotFound = .false.
       end if
    end do

    IF(NotFound) then
       IF (zlnew) THEN
          N_HANDLE     = N_HANDLE + 1

          IF(N_HANDLE .gt. NMAX_HANDLE) THEN
             CALL error_bi(&
                  'nmax_handle exceeded, recompile with higher nmax_handle'&
                  , substr)
          ENDIF

          ind              = N_HANDLE
          proc(ind)%name   = TRIM(ADJUSTL(name))
          proc(ind)%handle = ind

          DO id = 1, n_dom
             ALLOCATE(proc(ind)%p(id)%lvmod(NProgVmax:-1))
             DO jv = NProgVmax,-1
                proc(ind)%p(id)%lvmod(jv) = .FALSE.
             END DO
          END DO
       ELSE

          ind = -1

       END IF
    END IF

    MTEND_GET_HANDLE = ind

  END FUNCTION mtend_get_handle
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------
  SUBROUTINE mtend_register(handle, mtend_id, mtend_id_vec, dom_id)

    USE messy_main_channel_mem,  ONLY: dom_curid

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)                         :: handle   ! process handle
    INTEGER              , INTENT(IN), OPTIONAL :: mtend_id ! internal no. of prog. variables
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: mtend_id_vec ! tracer indices
    INTEGER              , INTENT(IN), OPTIONAL :: dom_id     !

    ! LOCAL
    CHARACTER(LEn=*), PARAMETER        :: substr  = 'mtend_register'
    INTEGER, DIMENSION(:), ALLOCATABLE :: idxtmp
    INTEGER                            :: n, m, pix

    IF (PRESENT(mtend_id) .AND. PRESENT(mtend_id_vec)) THEN
       CALL error_bi('call subroutine either with mtend_id_vec or mtend_id',substr)
    ENDIF

    IF (.NOT. PRESENT(mtend_id) .AND. .NOT. PRESENT(mtend_id_vec)) THEN
       CALL error_bi('subroutine call requires one of mtend_id_vec or mtend_id',substr)
    ENDIF

    IF (PRESENT(dom_id)) THEN
       pix = dom_id
    ELSE
       pix = dom_curid
    END IF

    ! NOTE: it is NOT tested, if tracer IDs occur more than once!
    IF (PRESENT(mtend_id)) THEN
       IF (mtend_id < 0) THEN
          proc(handle)%p(pix)%lvmod(mtend_id) = .TRUE.
          RETURN
       ELSE IF  (mtend_id == mtend_id_tracer) THEN
          m = 0
          IF (ASSOCIATED(proc(handle)%p(pix)%idx)) THEN
             CALL error_bi(&
                  'call subroutine only once for mtend_id_tracer '//TRIM(proc(handle)%name)&
                  ,substr)
          ENDIF
       ELSE
          m = 1
       END IF
    ELSE
       m = SIZE(mtend_id_vec)
    END IF

    IF (ASSOCIATED(proc(handle)%p(pix)%idx)) THEN
       n = SIZE(proc(handle)%p(pix)%idx)
       ALLOCATE(idxtmp(n))
       idxtmp(:) = proc(handle)%p(pix)%idx(:)
       DEALLOCATE(proc(handle)%p(pix)%idx) ; NULLIFY(proc(handle)%p(pix)%idx)
    ELSE
       n = 0
    END IF

    ALLOCATE(proc(handle)%p(pix)%idx(n+m))

    IF (ALLOCATED(idxtmp)) THEN
       proc(handle)%p(pix)%idx(1:n) = idxtmp(:)
       DEALLOCATE(idxtmp)
    ENDIF

    IF (PRESENT(mtend_id_vec)) THEN
       proc(handle)%p(pix)%idx(n+1:n+m) = mtend_id_vec(:)
    ELSE IF (m > 0) THEN
       proc(handle)%p(pix)%idx(n+1)     = mtend_id
    END IF

    proc(handle)%p(pix)%any_tracer = .TRUE.

  END SUBROUTINE mtend_register
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------
  SUBROUTINE mtend_get_start_l(mtend_id, v0, v0t)

    USE messy_main_grid_def_mem_bi,  ONLY: kproma,jrow

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER   :: substr = 'mtend_get_start_l'
    INTEGER, INTENT(IN)           :: mtend_id   ! variable identifier
    REAL(kind=dp), DIMENSION(:,:),   OPTIONAL, INTENT(OUT) :: v0 ! start value
    ! start value of all tracers
    REAL(kind=dp), DIMENSION(:,:,:), OPTIONAL, INTENT(OUT) :: v0t

    ! LOCAL
    INTEGER :: jt

    IF (mtend_id == mtend_id_tracer .AND. .NOT. PRESENT(v0t)) &
         CALL error_bi('marker mtend_id_tracer only correct with v0t', substr)

    IF (PRESENT(v0t)) THEN

       ! ALL_TRACERS
       v0t(:,:,:) = 0.0_dp
       v0t(_RI_X_ZN_(1:kproma,:,:)) = xtm1(_RI_XYZN_(1:kproma,jrow,:,:)) + &
            xtte(_RI_XYZN_(1:kproma,jrow,:,:))*time_step_len

    ELSE
       IF (progvar(mtend_id)%lexists) THEN

          v0(:,:) = 0.0_dp
          v0(1:kproma,:) = &
               progvar(mtend_id)%varm1ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
               progvar(mtend_id)%varteptr(_RI_XYZ__(1:kproma,jrow,:))*time_step_len
       ELSE
          write(iouerr,*) substr, 'requested variable:' , mtend_id
          CALL error_bi('requested variable does not exist, provide calculation routine'&
               , substr)
       END IF

    END IF !tracer

  END SUBROUTINE mtend_get_start_l
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------
  SUBROUTINE mtend_get_start_g(mtend_id, v0, v0t)

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER   :: substr = 'mtend_get_start_g'
    INTEGER, INTENT(IN)           :: mtend_id   ! variable identifier
    REAL(kind=dp), DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)   :: v0 ! start value
    ! start value of all tracers
    REAL(kind=dp), DIMENSION(:,:,:,:), INTENT(OUT), OPTIONAL :: v0t


    IF (mtend_id == mtend_id_tracer .AND. .NOT. PRESENT(v0t)) &
         CALL error_bi('marker mtend_id_tracer only correct with v0t', substr)

    IF (PRESENT(v0t)) THEN

       ! ALL_TRACERS
       v0t(:,:,:,:) = xtm1(:,:,:,:) +  xtte(:,:,:,:)*time_step_len

    ELSE
       IF (progvar(mtend_id)%lexists) THEN

          v0(:,:,:) = progvar(mtend_id)%varm1ptr(:,:,:) + &
               progvar(mtend_id)%varteptr(:,:,:)*time_step_len
       ELSE
          CALL error_bi('requested variable does not exist, provide calculation routine'&
               , substr)
       END IF

    END IF

  END SUBROUTINE mtend_get_start_g
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------
  SUBROUTINE mtend_add_l (handle, mtend_id, px, pxt, l_add)

    USE messy_main_grid_def_mem_bi,    ONLY: jrow, kproma, nproma, nlev
    USE messy_main_tools,              ONLY: int2str ! op_pj_20160801
    USE messy_main_channel_mem,        ONLY: dom_curid

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*), PARAMETER :: substr = 'mtend_add_l'
    INTEGER, INTENT(IN) :: handle    ! process identifier
    INTEGER, INTENT(IN) :: mtend_id  ! variable identifier
    REAL(kind=dp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: px ! tendency to add
    ! tendency of all tracers
    REAL(kind=dp), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: pxt
    LOGICAL, INTENT(IN), OPTIONAL                         :: l_add

    ! LOCAL
    INTEGER          :: jt, jk, ic, pix
    REAL(kind=dp), DIMENSION(:,:), POINTER :: zz => NULL()
    CHARACTER(LEN=5) :: idtstr = '     '
    LOGICAL          :: ladd     ! add to tendency field

    !-------------------------
    !for development diagnostics
    !-------------------------
    REAL(kind=dp), DIMENSION(:,:,:), POINTER :: z3d => NULL()
    REAL(kind=dp), DIMENSION(:,:),   POINTER :: z2d => NULL()
    !-------------------------

    !write(iouerr,*) substr, ' ',TRIM(proc(handle)%name),dom_curid,' ', mtend_id

    ! SET LEVEL DIMENSION
    IF (PRESENT(px) .AND. PRESENT(pxt)) THEN
       CALL error_bi('only one of px or pxt must be present for call',substr)
    END IF
    IF (mtend_id == mtend_id_tracer .AND. .NOT. PRESENT(pxt)) &
         CALL error_bi('marker mtend_id_tracer only correct with pxt', substr)

    IF (PRESENT(l_add)) THEN
       ladd = l_add
    ELSE
       ladd = .TRUE.
    ENDIF

    IF (mtend_id /= mtend_id_w) THEN
       ALLOCATE(zz(nproma,nlev))
    ELSE
       ALLOCATE(zz(nproma,nlev+1))
    END IF

    ! Internal horic. wind vector tendency is never scaled with cos lat.
    ! This is to create factor for correct scaling of incoming tendency
    IF((wind_is_scaled) .AND. ((mtend_id == mtend_id_u) .OR. &
         (mtend_id == mtend_id_v)))  THEN
       DO jk=1,SIZE(zz,2)
#ifdef WINDSCALING
          zz(:,jk) = sqcst_2d(:,jrow)
#else
          zz(:,jk) = 1._dp
#endif
          ! Set the zeros of the remaining piece (nproma-npromz)
          ! to 1 to avoid getting floating point error
          where (zz(:,jk) < 1.0E-40_dp)
             zz(:,jk) = 1._dp
          end where
       END DO
    ELSE
       zz(:,:) = 1._dp
    ENDIF

    !-------------------------

    IF (mtend_id < NProgVmax .OR. mtend_id > proc(handle)%p(dom_curid)%ntrac_proc) THEN
       write(iouerr,*) 'ERROR in TENDENCY', mtend_id, '> ' &
            , proc(handle)%p(dom_curid)%ntrac_proc, ' .OR. < ', NProgVmax
       write(iouerr,*) 'ERROR in TENDENCY: handle', handle, TRIM(proc(handle)%name)
       CALL error_bi(' VARIABLE IDENTIFIER OUT OF RANGE',substr)
     END IF

    ! ALL TRACERS
    IF (PRESENT(pxt)) THEN
       IF (proc(handle)%p(dom_curid)%ntrac_proc < 0) &
            CALL error_bi(' TRACERS UNREGISTERED IN '//&
            &TRIM(proc(handle)%name), substr)

       tracer_loop: DO jt=1, proc(handle)%p(dom_curid)%ntrac_proc
          IF (.NOT. proc(handle)%p(dom_curid)%lvmod(jt)) THEN
             CALL int2str(idtstr,mtend_id,' ','X')
             CALL error_bi(' UNREGISTERED TRACER '//TRIM(idtstr)//' IN '//&
                  &TRIM(proc(handle)%name),substr)
          END IF

          ! UPDATE "total internal" tendency
          IF (ladd) THEN
             xtte(_RI_XYZN_(1:kproma,jrow,:,jt)) =    &
                  xtte(_RI_XYZN_(1:kproma,jrow,:,jt)) &
                  + pxt( _RI_X_ZN_(1:kproma,:,jt))
          END IF

          ! UPDATE sum of tendencies
          IF(l_closure(dom_curid)) THEN
             IF (handle /= I_HANDLE_BMTT) THEN
                proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
                     proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:))&
                     + pxt( _RI_X_ZN_(1:kproma,:,jt))
             END IF
          ENDIF

          ! SAVE COPY OF TENDENCY IN CHANNEL OBJECTS
          DO ic = 1, NCASES
             IF (proc(handle)%p(dom_curid)%case(ic)%lvar(jt)) THEN
                proc(handle)%p(dom_curid)%case(ic)%pvar(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
                     proc(handle)%p(dom_curid)%case(ic)%pvar(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:)) &
                     + pxt(_RI_X_ZN_(1:kproma,:,jt))
             ENDIF
          ENDDO

       END DO tracer_loop
       !for development diagnostics
       IF (l_clos_diag(dom_curid)) z3d => xtte(_RI_XYZN_(1:kproma,jrow,:,:))

    ELSE
       ! "ordinary" PROGNOSTIC VARIABLE or single TRACERS
       IF (mtend_id > 0 .AND. .NOT. proc(handle)%p(dom_curid)%ntrac_proc > 0) &
            CALL error_bi(' TRACERS UNREGISTERED IN '//&
            &TRIM(proc(handle)%name), substr)
       IF (.NOT. proc(handle)%p(dom_curid)%lvmod(mtend_id)) THEN
          CALL int2str(idtstr,mtend_id,' ','X')
          CALL error_bi(' UNREGISTERED TENDENCY '//TRIM(idtstr)//' IN '//&
               &TRIM(proc(handle)%name),substr)
       ENDIF

        ! UPDATE "total internal" tendency
       IF (ladd) THEN
          progvar(mtend_id)%varteptr(_RI_XYZ__(1:kproma,jrow,:)) = &
               progvar(mtend_id)%varteptr(_RI_XYZ__(1:kproma,jrow,:)) + px(1:kproma,:)
          IF (l_clos_diag(dom_curid))  &
               z2d => progvar(mtend_id)%varteptr(_RI_XYZ__(1:kproma,jrow,:))
       END IF

       ! UPDATE sum of tendencies
       IF(l_closure(dom_curid)) THEN
          IF (handle /= I_HANDLE_BMTT) THEN ! ub_ak_20190716
             proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(mtend_id)%ptr(_RI_XYZ__(1:kproma,jrow,:)) =    &
                  proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(mtend_id)%ptr(_RI_XYZ__(1:kproma,jrow,:)) &
                  + (px(1:kproma,:)/zz(1:kproma,:))
          END IF
       ENDIF !remove scaling if incoming tendency is scaled

       ! SAVE COPY OF TENDENCY IN CHANNEL OBJECTS
       DO ic = 1, NCASES
          IF (proc(handle)%p(dom_curid)%case(ic)%lvar(mtend_id)) THEN
             proc(handle)%p(dom_curid)%case(ic)%pvar(mtend_id)%ptr(_RI_XYZ__(1:kproma,jrow,:)) =    &
                  proc(handle)%p(dom_curid)%case(ic)%pvar(mtend_id)%ptr(_RI_XYZ__(1:kproma,jrow,:)) &
                  + (px(1:kproma,:)/zz(1:kproma,:))
          END IF
       END DO

    END IF

 !-------------------------
 !for development diagnostics
 !-------------------------
 IF (l_closure(dom_curid) .AND. l_clos_diag(dom_curid) .AND. ladd) THEN

       IF (PRESENT(pxt)) THEN
       DO jt =1, proc(handle)%p(dom_curid)%ntrac_proc
          write(*,*) p_pe, '=pe,',' ',nump,' ','LOCAL_ADD',' ','All_Tracers!.. Process handle:', handle &
               , TRIM(proc(handle)%name),' ','Tracer Nr:',jt &
               ,'min(add_tend)=',minval(pxt(_RI_X_ZN_(1:kproma,:,jt))),'max(add_tend)=',maxval(pxt(_RI_X_ZN_(1:kproma,:,jt))) &
               ,'min(ext_sum_tend)=',minval(z3d(_RI_X_ZN_(1:kproma,:,jt))),'max(ext_sum_tend)=',maxval(z3d(_RI_X_ZN_(1:kproma,:,jt))) &
               ,'min(int_sum_tend)=',MINVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:)))&
               ,'max(int_sum_tend)=',MAXVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:))) &
               ,'min(int-ext)=',MINVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:))-z3d(_RI_X_ZN_(1:kproma,:,jt))) &
               ,'max(int-ext)=',MAXVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:))-z3d(_RI_X_ZN_(1:kproma,:,jt)))
       ENDDO
    ELSE
       write(*,*) p_pe, '=pe,',' ',nump,' ','LOCAL_ADD',' ','Prog_Variable!.. Process handle:', handle &
            , TRIM(proc(handle)%name),' ','Variable:',mtend_id &
            ,'min(add_tend)=',minval(px(1:kproma,:)/zz(1:kproma,:)),'max(add_tend)=',maxval(px(1:kproma,:)/zz(1:kproma,:)) &
            ,'min(ext_sum_tend)=',minval(z2d(1:kproma,:)/zz(1:kproma,:))&
            ,'max(ext_sum_tend)=',maxval(z2d(1:kproma,:)/zz(1:kproma,:))&
            ,'min(int_sum_tend)=',MINVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(mtend_id)%ptr(_RI_XYZ__(1:kproma,jrow,:)))&
            ,'max(int_sum_tend)=',MAXVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(mtend_id)%ptr(_RI_XYZ__(1:kproma,jrow,:)))&
            ,'min(int-ext)=',MINVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(mtend_id)%ptr(_RI_XYZ__(1:kproma,jrow,:))-z2d(1:kproma,:)/zz(1:kproma,:))&
            ,'max(int-ext)=',MAXVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(mtend_id)%ptr(_RI_XYZ__(1:kproma,jrow,:))-z2d(1:kproma,:)/zz(1:kproma,:))
    ENDIF
    write(*,*)'- - - - - - - - - - - - - - -'
    nump=nump+1
 ENDIF

 DEALLOCATE(zz)
 NULLIFY(zz)
 !-------------------------

  END SUBROUTINE mtend_add_l
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------
  SUBROUTINE mtend_add_g (handle, mtend_id, px, pxt, l_add)

    USE messy_main_grid_def_mem_bi,    ONLY: nproma, ngpblks, nlev
    USE messy_main_tools,              ONLY: int2str
    USE messy_main_channel_mem,        ONLY: dom_curid

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER                           :: substr = 'mtend_add_g'
    INTEGER, INTENT(IN)                                   :: handle   ! process identifier
    INTEGER, INTENT(IN)                                   :: mtend_id ! variable identifier
    REAL(kind=dp), DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: px       ! tendency to add

    ! tendency of all tracers
    REAL(kind=dp), DIMENSION(:,:,:,:), INTENT(IN), OPTIONAL                :: pxt
    LOGICAL, INTENT(IN), OPTIONAL                         :: l_add

    ! LOCAL
    INTEGER :: jt, jk, ic
    CHARACTER(LEN=5) :: idtstr = '     '
    REAL(kind=dp), DIMENSION(:,:,:), POINTER :: zz => NULL()
    LOGICAL :: ladd     ! add to tendency field

    !-------------------------
    !for development diagnostics
    !-------------------------
    REAL(kind=dp), DIMENSION(:,:,:,:), POINTER :: z4d => NULL()
    REAL(kind=dp), DIMENSION(:,:,:), POINTER   :: z3d => NULL()
    !-------------------------

    !write(iouerr,*) substr, ' ',TRIM(proc(handle)%name), dom_curid,' ', mtend_id

    ! SET LEVEL DIMENSION
    IF (PRESENT(px) .AND. PRESENT(pxt)) THEN
       CALL error_bi('only one of px or pxt must be present for call',substr)
    END IF

    IF (mtend_id == mtend_id_tracer .AND. .NOT. PRESENT(pxt)) &
         CALL error_bi('marker mtend_id_tracer only correct with pxt', substr)

    IF (mtend_id < NProgVmax - 1 .OR. &
         mtend_id > proc(handle)%p(dom_curid)%ntrac_proc) THEN
       write(iouerr,*) 'ERROR in TENDENCY', mtend_id, '> ' &
            , proc(handle)%p(dom_curid)%ntrac_proc, ' .OR. < ', NProgVmax
       write(iouerr,*) 'ERROR in TENDENCY: handle', handle, TRIM(proc(handle)%name)
       CALL error_bi(' VARIABLE IDENTIFIER OUT OF RANGE',substr)
    END IF

    IF (PRESENT(l_add)) THEN
       ladd = l_add
    ELSE
       ladd = .TRUE.
    ENDIF

    IF (mtend_id /= mtend_id_w) THEN
       ALLOCATE(zz(_RI_XYZ__(nproma,ngpblks,nlev)))
    ELSE
       ALLOCATE(zz(_RI_XYZ__(nproma,ngpblks,nlev+1)))
    END IF

    ! Internal horic. wind vector tendency is never scaled with cos lat.
    ! This is to create factor for correct scaling of incoming tendency
    IF((wind_is_scaled) .AND. ((mtend_id == mtend_id_u) .OR. (mtend_id == mtend_id_v)))  THEN
       DO jk=1,SIZE(zz, _IZ_XYZ__ )
#ifdef WINDSCALING
          zz(_RI_XYZ__(:,:,jk)) = sqcst_2d(:,:)
#else
          zz(_RI_XYZ__(:,:,jk)) = 1._dp
#endif
          ! Set the zeros of the remaining piece (nproma-npromz)
          ! to 1 to avoid getting floating point error
          where(zz(_RI_XYZ__(:,:,jk)) < 1.0E-40)
             zz(_RI_XYZ__(:,:,jk)) = 1._dp
          end where
        END DO
    ELSE
       zz(:,:,:) = 1._dp
    ENDIF
    !-------------------------

    ! ALL TRACERS
    IF (PRESENT(pxt)) THEN
       IF (proc(handle)%p(dom_curid)%ntrac_proc < 0) &
            CALL error_bi(' TRACERS UNREGISTERED IN '//&
            &TRIM(proc(handle)%name), substr)

       tracer_loop: DO jt=1,  proc(handle)%p(dom_curid)%ntrac_proc
          IF (.NOT. proc(handle)%p(dom_curid)%lvmod(jt)) THEN
             CALL int2str(idtstr,mtend_id,' ','X')
             CALL error_bi(' UNREGISTERED TRACER '//TRIM(idtstr)//' IN '//&
                  &TRIM(proc(handle)%name),substr)
          END IF
          ! UPDATE "total internal" tendency
          IF (ladd) &
               xtte(_RI_XYZN_(:,:,:,jt)) = xtte(_RI_XYZN_(:,:,:,jt)) + pxt(_RI_XYZN_(:,:,:,jt))
          ! UPDATE sum of tendencies
          IF(l_closure(dom_curid)) THEN
             IF (handle /= I_HANDLE_BMTT) THEN
                proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jt)%ptr(:,:,:) = &
                     proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jt)%ptr(:,:,:) &
                     + pxt(_RI_XYZN_(:,:,:,jt))
             END IF
          ENDIF
          ! SAVE COPY OF TENDENCY IN CHANNEL OBJECTS
          DO ic = 1, NCASES
             IF (proc(handle)%p(dom_curid)%case(ic)%lvar(jt)) THEN
                proc(handle)%p(dom_curid)%case(ic)%pvar(jt)%ptr(:,:,:) = &
                     proc(handle)%p(dom_curid)%case(ic)%pvar(jt)%ptr(:,:,:)&
                     + pxt(_RI_XYZN_(:,:,:,jt))
             END IF
          END DO
       ENDDO tracer_loop

       IF (l_clos_diag(dom_curid)) z4d => xtte   !for development diagnostics(:,:,jt,:)

    ELSE

       ! "ordinary" PROGNOSTIC VARIABLE or single TRACERS
       IF (mtend_id > 0 .AND. .NOT. proc(handle)%p(dom_curid)%any_tracer) &
            CALL error_bi(' TRACERS UNREGISTERED IN '//&
            &TRIM(proc(handle)%name), substr)
       IF (.NOT. proc(handle)%p(dom_curid)%lvmod(mtend_id)) THEN
          CALL int2str(idtstr,mtend_id,' ','X')
          CALL error_bi(' UNREGISTERED TENDENCY '//TRIM(idtstr)//' IN '//&
                  &TRIM(proc(handle)%name),substr)
       END IF

       ! UPDATE "total global" tendency
       IF (ladd) THEN
          progvar(mtend_id)%varteptr = progvar(mtend_id)%varteptr + px(:,:,:)
          IF (l_clos_diag(dom_curid)) z3d => progvar(mtend_id)%varteptr
       END IF

       ! UPDATE sum of tendencies
       IF(l_closure(dom_curid)) THEN
          IF (handle /= I_HANDLE_BMTT) THEN
             proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(mtend_id)%ptr(:,:,:) = &
                  proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(mtend_id)%ptr(:,:,:) +&
                  px(:,:,:)/zz(:,:,:)
          END IF
       ENDIF

       ! SAVE COPY OF TENDENCY IN CHANNEL OBJECTS
       DO ic = 1, NCASES
          IF (proc(handle)%p(dom_curid)%case(ic)%lvar(mtend_id)) THEN
             proc(handle)%p(dom_curid)%case(ic)%pvar(mtend_id)%ptr(:,:,:) = &
                  proc(handle)%p(dom_curid)%case(ic)%pvar(mtend_id)%ptr(:,:,:) + px(:,:,:)/zz(:,:,:)
          ENDIF
       END DO
    END IF

    !-----------------------
    !for development diagnostics
    IF (l_closure(dom_curid) .AND. l_clos_diag(dom_curid) .AND. ladd) THEN
       IF (PRESENT(pxt)) THEN
          DO jt = 1 ,  proc(handle)%p(dom_curid)%ntrac_proc
             write(*,*) p_pe, '=pe,',' ',nump,' ','GLOBAL_ADD',' ','All_Tracers!.. Process handle:', handle &
                  ,TRIM(proc(handle)%name),' ','Tracer Nr:',jt &
                  ,'min(add_tend)=' , minval(pxt(_RI_XYZN_(:,:,:,jt))),'max(add_tend)=', maxval(pxt(_RI_XYZN_(:,:,:,jt))) &
                  ,'min(ext_sum_tend)=', minval(z4d(_RI_XYZN_(:,:,:,jt))),'max(ext_sum_tend)=', maxval(z4d(_RI_XYZN_(:,:,:,jt))) &
                  ,'min(int_sum_tend)=',MINVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jt)%ptr)&
                  ,'max(int_sum_tend=)',MAXVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jt)%ptr) &
                  ,'min(int-ext)=',MINVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jt)%ptr-z4d(_RI_XYZN_(:,:,:,jt)))&
                  ,'max(int-ext)=',MAXVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(jt)%ptr-z4d(_RI_XYZN_(:,:,:,jt)))
          ENDDO
       ELSE
          write(*,*) p_pe, '=pe,',' ',nump,' ','GLOBAL_ADD',' ','Prog_Variable!.. Process handle:',handle &
               , TRIM(proc(handle)%name),' ','Variable:',mtend_id &
               ,'min(add_tend)=',minval(px(:,:,:)/zz(:,:,:)),'max(add_tend)=',maxval(px(:,:,:)/zz(:,:,:)) &
               ,'min(ext_sum_tend)=', minval(z3d/zz(:,:,:)),'max(ext_sum_tend)=', maxval(z3d/zz(:,:,:)) &
               ,'min(int_sum_tend)=',MINVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(mtend_id)%ptr)&
               ,'max(int_sum_tend)=',MAXVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(mtend_id)%ptr)&
               ,'min(int-ext)=',MINVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(mtend_id)%ptr-z3d/zz(:,:,:))&
               ,'max(int-ext)=',MAXVAL(proc(I_HANDLE_SUM)%p(dom_curid)%case(IC_GEN)%pvar(mtend_id)%ptr-z3d/zz(:,:,:))
       ENDIF
       write(*,*)'- - - - - - - - - - - - - - -'
       nump=nump+1
    ENDIF

    DEALLOCATE(zz)
    NULLIFY(zz)
    !-----------------------

  END SUBROUTINE mtend_add_g
  ! ------------------------------------------------------------------

  !-------------------------------------------------------------------
  SUBROUTINE mtend_request(process_string, mtend_id, ptr_out, dom_id &
       , lextern)

    USE messy_main_channel,          ONLY: new_channel_object, new_attribute
    USE messy_main_channel_mem,      ONLY: dom_curid
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_tracer,           ONLY: get_tracer, STRLEN_TRSET
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    USE messy_main_constants_mem,    ONLY: STRLEN_MEDIUM
#ifdef ICON
    USE messy_main_tracer_mem_bi,    ONLY: L_GPTRSTR
#endif

    IMPLICIT NONE

    CHARACTER(len=*),  INTENT(IN)                    :: process_string
    INTEGER,           INTENT(IN)                    :: mtend_id
    REAL(KIND=DP),     POINTER,    DIMENSION(:,:,:)  :: ptr_out
    INTEGER,           INTENT(IN), OPTIONAL          :: dom_id
    LOGICAL,           INTENT(IN), OPTIONAL          :: lextern

    CHARACTER(LEN=*), PARAMETER               :: substr = 'mtend_request'
    INTEGER                                   :: status
    INTEGER                                   :: exch_handle
    CHARACTER(LEN=2*STRLEN_MEDIUM+1)          :: name     = ''
    CHARACTER(LEN=2*STRLEN_MEDIUM+1)          :: name1    = ''
    CHARACTER(LEN=STRLEN_MEDIUM)              :: unit     = ''
    CHARACTER(LEN=STRLEN_MEDIUM)              :: unit1    = ''
    INTEGER                                   :: pix
    CHARACTER(LEN=STRLEN_TRSET)               :: GP_TRSTR
    LOGICAL                                   :: l_extern

    IF (PRESENT(dom_id)) THEN
       pix = dom_id
    ELSE
       pix = dom_curid
    END IF
    IF (PRESENT(lextern)) THEN
       l_extern = lextern
    ELSE
       l_extern = .FALSE.
    END IF

    ! get requested process handle for exchange
    exch_handle = mtend_get_handle(process_string, .FALSE.)

    ! check if process is available
    IF (exch_handle .lt. 0) THEN
       CALL error_bi( &
            'process not existent, see process string in exchanging module: '//TRIM(process_string), substr)
    ENDIF

    IF (mtend_id == 0) CALL error_bi( ' 0  is an undefined variable identifier', substr)

    IF (l_extern) THEN
       IF ( proc(exch_handle)%p(pix)%case(IC_GEN)%lextern(mtend_id)) &
            CALL error_bi('MEMORY REQUESTED FOR EXTERN IS ALREADY EXTERN' &
            , substr)
    END IF

    TRACER:IF (mtend_id > mtend_id_tracer) THEN

#ifndef ICON
       GP_TRSTR = GPTRSTR
#else
       GP_TRSTR = L_GPTRSTR(pix)
#endif
       CALL get_tracer(status, GP_TRSTR, mtend_id, fullname = name1, unit = unit1)
       CALL tracer_halt(substr, status)
       name = TRIM(name1)
       unit = '('//TRIM(unit1)//') / s'
    ELSE
       name = pvarname(mtend_id)
       unit = pvarunit(mtend_id)
    END IF TRACER

    IF(.NOT. proc(exch_handle)%p(pix)%case(IC_GEN)%lvar(mtend_id))THEN
       !:::::::::::::::::::::::::::::::::::::::::::::::::::::
       ! create channel object for requested combination of process and tracer
       !:::::::::::::::::::::::::::::::::::::::::::::::::::::
       CALL new_channel_object(status, modstr,                     &
            TRIM(process_string)//'_'//TRIM(name),                      &
            p3 = proc(exch_handle)%p(pix)%case(IC_GEN)%pvar(mtend_id)%ptr, &
            reprid=progvar(mtend_id)%reprid                        &
#ifdef ICON
            , dom_id=pix                                           &
#endif
            )
       CALL channel_halt(substr, status)

       CALL new_attribute(status, modstr,                                      &
            TRIM(process_string)//'_'//TRIM(name),                                  &
            'long_name', c='Exchange Process'//TRIM(process_string)//TRIM(name)&
#ifdef ICON
            , dom_id=pix                            &
#endif
            )
       CALL channel_halt(substr, status)

       CALL new_attribute(status, modstr,                   &
            TRIM(process_string)//'_'//TRIM(name),   &
            'units', c='('//TRIM(unit)//') / s'             &
#ifdef ICON
            , dom_id=pix                &
#endif
            )
       CALL channel_halt(substr, status)
       !:::::::::::::::::::::::::::::::::::::::::::::::::::::
       proc(exch_handle)%p(pix)%case(IC_GEN)%lvar(mtend_id) = .TRUE.
    ENDIF

    ptr_out => proc(exch_handle)%p(pix)%case(IC_GEN)%pvar(mtend_id)%ptr(:,:,:)

    ! SET INDICATOR THAT FIELD CONTENT IS MANAGED EXTERN
    IF (l_extern) THEN
       proc(exch_handle)%p(pix)%case(IC_GEN)%lextern(mtend_id) = .TRUE.
       proc(exch_handle)%p(pix)%case(IC_GEN)%lvar(mtend_id)    = .FALSE.
    END IF

  END SUBROUTINE mtend_request
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------
  SUBROUTINE mtend_set_sqcst_scal (scale_flag)

    implicit    none
    logical,intent(IN)                     :: scale_flag

    ! Set Module in sqcst scaling mode
    ! This is done at begin and switched back to
    !  non scaling mode at end of physc

    wind_is_scaled = scale_flag

    RETURN

  END SUBROUTINE mtend_set_sqcst_scal
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------
  LOGICAL FUNCTION mtend_checkreg (process_string, mtend_id, dom_id)

    ! MESSy
    USE messy_main_channel_mem,       ONLY: dom_curid

    IMPLICIT NONE

    ! I/O
    CHARACTER(len=*),  INTENT(IN)                    :: process_string
    INTEGER,           INTENT(IN)                    :: mtend_id
    INTEGER,           INTENT(IN), OPTIONAL          :: dom_id

    ! LOCAL
    LOGICAL                                          :: flag
    INTEGER                                          :: handle = 0
    CHARACTER(LEN=*), PARAMETER                      :: substr = 'mtend_checkreg'
    INTEGER                                          :: pix

    IF (PRESENT(dom_id)) THEN
       pix = dom_id
    ELSE
       pix = dom_curid
    END IF

    ! Flag by default false
    flag = .FALSE.

    ! get requested process handle
    handle = mtend_get_handle(process_string, .FALSE.)

    IF (handle .ge. 0) THEN
       TRACER:IF (mtend_id > mtend_id_tracer) THEN

          IF (proc(handle)%p(pix)%any_tracer) THEN
             IF (SIZE(proc(handle)%p(pix)%idx) .le. 0) THEN
                flag = .TRUE.
             ENDIF

             IF (SIZE(proc(handle)%p(pix)%idx) .gt. 0) THEN
                IF (ANY(proc(handle)%p(pix)%idx .eq. mtend_id)) THEN
                   flag = .TRUE.
                ENDIF
             ENDIF
          ENDIF

       ELSE

          flag = proc(handle)%p(pix)%lvmod(mtend_id)

       ENDIF TRACER
    ENDIF

    MTEND_CHECKREG = flag

  END FUNCTION mtend_checkreg
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------
  SUBROUTINE tendency_read_nml_cpl(status, iou)


    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ l_closure, l_clos_diag, l_full_diag, TDIAG

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='tendency_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

    return

  END SUBROUTINE tendency_read_nml_cpl
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  SUBROUTINE tendency_parse_nml_cpl

    USE messy_main_tools,            ONLY: strcrack, int2str
    USE messy_main_channel,          ONLY: new_channel_object, new_attribute
    USE messy_main_channel_error_bi, ONLY: channel_halt
#ifdef ICON
    USE messy_main_channel_bi,       ONLY: main_channel_set_domain
    USE messy_main_tracer_tools_bi,  ONLY: main_tracer_set_domain
#endif
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    USE messy_main_constants_mem,    ONLY: STRLEN_MEDIUM
    USE messy_main_tracer,           ONLY: get_tracer,STRLEN_TRSET

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER       :: substr='tendency_parse_nml_cpl'

    INTEGER                           :: jd, nr1, nr2, ih, jc, jv
    CHARACTER(len=132),DIMENSION(:),POINTER   :: Pname1 => NULL()
    CHARACTER(len=132),DIMENSION(:),POINTER   :: Pname2 => NULL()

    CHARACTER(len=3)                  :: int_string
    CHARACTER(len=STRLEN_MEDIUM)      :: sumunit
    CHARACTER(LEN=2*STRLEN_MEDIUM+1)  :: tracname = ''
    CHARACTER(LEN=STRLEN_MEDIUM)      :: basename = ''
    CHARACTER(LEN=STRLEN_MEDIUM)      :: subname  = ''
    CHARACTER(LEN=STRLEN_MEDIUM)      :: unit     = ''
    INTEGER                           :: si
    CHARACTER(LEN=STRLEN_TRSET)       :: GP_TRSTR = ''

    LOGICAL :: NoVariable

    !loop for each variable in namelist
    tdiag_loop: DO jd=1,N_TDIAG

       !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
       !To assign unit of variable to channel object
       !``````````````````````````````````````````````````
       NoVariable = .TRUE.

#ifdef ICON
       CALL main_channel_set_domain(XTDIAG(jd)%io%domain_idx)
       CALL main_tracer_set_domain(XTDIAG(jd)%io%domain_idx)
       CALL main_tendency_set_domain(XTDIAG(jd)%io%domain_idx)
#endif

       var_loop: DO  jv = NProgVmax, -1

          IF(TRIM(ADJUSTL(XTDIAG(jd)%io%variable_name))&
               .eq. TRIM(ADJUSTL(pvarname(jv)))) THEN
             sumunit = pvarunit(jv)
             NoVariable = .FALSE.
             XTDIAG(jd)%progvaridx = jv
             EXIT ! use jv for channel object definition
          ENDIF

       ENDDO var_loop

       ! be careful: tracers are not in the list of variables ...
       IF (NoVariable) THEN
          ! if variable is a tracer:
          tracname=ADJUSTL(TRIM(XTDIAG(jd)%io%variable_name))
          si = INDEX(tracname, '_')
          IF (si == 0) THEN
             basename = TRIM(tracname)
             subname  = ''
          ELSE
             basename = TRIM(tracname(:si-1))
             subname  = TRIM(tracname(si+1:))
          END IF
          CALL get_tracer(status, GPTRSTR, basename, subname=subname &
               , unit=unit, idx=XTDIAG(jd)%progvaridx)
          IF (status == 0) THEN
             ! TRACER EXISTS
             sumunit='('//TRIM(ADJUSTL(unit))//')/s'
          ELSE
             ! TRACER DOES NOT EXIST
             CALL error_bi(&
                  'requested variable or tracer in tendency nml not existent: '&
                  //TRIM(ADJUSTL(XTDIAG(jd)%io%variable_name)), substr)
          ENDIF
       END IF
       !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
       !``````````````````````````````````````````````````

       ! first crack whole process-string into seperate sums
       call strcrack (XTDIAG(jd)%io%process_name, ';', Pname1, nr1)

       ! number of sums (classes) of processes plus one for the rest
       XTDIAG(jd)%ntclass = nr1 + 1

       ! nproc = number of processes per sum
       ALLOCATE(XTDIAG(jd)%nproc(XTDIAG(jd)%ntclass))
       ! allocate pointer array for channel objects
       ALLOCATE(XTDIAG(jd)%obj(XTDIAG(jd)%ntclass))
       ! space for handles of processes
       ALLOCATE(XTDIAG(jd)%iproc(XTDIAG(jd)%ntclass))

       ! loop over sums for each nml entry plus one for rest
       class_loop: DO jc = 1, XTDIAG(jd)%ntclass

          !convert integer 'jc' into string
          CALL int2str(int_string,jc)

          !:::::::::::::::::::::::::::::::::::::::::::::::::::::

          IF (jc < XTDIAG(jd)%ntclass) THEN
             ! crack sums of processes into seperate processes (handles)
             call strcrack (Pname1(jc), '+', Pname2, nr2)!Pname1(jd)

             !allocate iproc with number of processes per sum (class)
             ALLOCATE(XTDIAG(jd)%iproc(jc)%ptr(nr2))
             XTDIAG(jd)%nproc(jc) = nr2
             XTDIAG(jd)%iproc(jc)%ptr(:) = 0

             !loop to assign process handle to iproc
             handle_loop: DO ih = 1, nr2

                !assign handle of processes to iproc%ptr
                XTDIAG(jd)%iproc(jc)%ptr(ih) = mtend_get_handle(Pname2(ih),.FALSE.)

                IF(XTDIAG(jd)%iproc(jc)%ptr(ih) .lt. 0)THEN
                   CALL error_bi( &
                        'requested process in tendency nml not existent: '//TRIM(Pname2(ih)), substr)
                ENDIF

             ENDDO handle_loop

             if (ASSOCIATED(Pname2)) deallocate(Pname2)

          ELSE
             ! UNACCOUNTED
             ALLOCATE(XTDIAG(jd)%iproc(jc)%ptr(1))
             XTDIAG(jd)%nproc(jc) = 1
             XTDIAG(jd)%iproc(jc)%ptr(:) = 0
          END IF

          !:::::::::::::::::::::::::::::::::::::::::::::::::::::
          ! create channel object for each sum (to be pointed on by process channels)
          !:::::::::::::::::::::::::::::::::::::::::::::::::::::

          CALL new_channel_object(status, modstr//chpf(IR_NML), &
               TRIM(XTDIAG(jd)%io%variable_name)//'_nml_sum_nr'//TRIM(int_string), &
               p3 = XTDIAG(jd)%obj(jc)%ptr,                                    &
               reprid=progvar(XTDIAG(jd)%progvaridx)%reprid                    &
               )
          CALL channel_halt(substr, status)

          IF (jc < XTDIAG(jd)%ntclass) THEN
             IF (XTDIAG(jd)%nproc(jc) > 1) THEN
                CALL new_attribute(status,modstr//chpf(IR_NML),&
                     TRIM(XTDIAG(jd)%io%variable_name)//'_nml_sum_nr'//TRIM(int_string),&
                     'long_name', c=TRIM(XTDIAG(jd)%io%variable_name)//': Sum of tendencies: '//TRIM(Pname1(jc))&
#ifdef ICON
                     , dom_id=XTDIAG(jd)%io%domain_idx &
#endif
                     )
             ELSE
                CALL new_attribute(status,modstr//chpf(IR_NML),&
                     TRIM(XTDIAG(jd)%io%variable_name)//'_nml_sum_nr'//TRIM(int_string),&
                     'long_name', c=TRIM(XTDIAG(jd)%io%variable_name)//': tendencies: '//TRIM(Pname1(jc))&
#ifdef ICON
                     , dom_id=XTDIAG(jd)%io%domain_idx &
#endif
                     )
             END IF
          ELSE
             CALL new_attribute(status,modstr//chpf(IR_NML),&
                  TRIM(XTDIAG(jd)%io%variable_name)//'_nml_sum_nr'//TRIM(int_string),&
                  'long_name', c='Sum of unaccounted tendencies'&
               )
          END IF
          CALL channel_halt(substr, status)


          CALL new_attribute(status,modstr//chpf(IR_NML),&
               TRIM(XTDIAG(jd)%io%variable_name)//'_nml_sum_nr'//TRIM(int_string),&
               'units', c=TRIM(sumunit)                  &
               )
          CALL channel_halt(substr, status)

       ENDDO class_loop

       if (ASSOCIATED(Pname1)) deallocate(Pname1)

    ENDDO tdiag_loop

    RETURN
  END SUBROUTINE tendency_parse_nml_cpl
  !-------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE compute_eps_and_clear (text,array,array_diff)

#ifdef COSMO
    USE messy_main_grid_def_mem_bi, ONLY: istart,iend, jstart, jend, ie, je
#endif

    IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)                          :: text
    REAL(kind=dp),DIMENSION(:,:,:),INTENT(IN)            :: array
    REAL(kind=dp),DIMENSION(:,:,:),INTENT(INOUT)         :: array_diff
    REAL(kind=dp)                                        :: eps,error
    INTEGER                                              :: ix

    eps = maxval(abs(array)) * 1.0E-10

    where(abs(array_diff) < eps)
       array_diff = 0.0
    end where
#ifdef COSMO
    array_diff(1:istart-1,:,:)   = 0.0
    array_diff(iend+1:ie,:,:)    = 0.0
    array_diff(:,1:jstart-1,:)   = 0.0
    array_diff(:,jend+1:je,:)    = 0.0
#endif
    error = maxval(abs(array_diff))

    if(error > eps)   then
       write(*,*) '==> WARNING: Difference in tendency too large:  ', &
            trim(text),p_pe,error,eps,MINVAL(array),MAXVAL(array),MINVAL(array_diff),MAXVAL(array_diff)
    end if

    return

  END SUBROUTINE compute_eps_and_clear
  ! ---------------------------------------------------------------

  !###########################################################################
  !---------------------------------------------------------------------------
  !###########################################################################
  ! ub_ak_20180928+
  SUBROUTINE set_basemodels_progvarset

    ! MESSy/BNIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! MESSy/SMCL
    USE messy_main_channel,          ONLY: get_channel_object
    USE messy_main_channel_repr,     ONLY: get_representation_id
#ifdef COSMO
    USE messy_main_channel,          ONLY: new_channel_object &
                                         , new_attribute
#endif

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'tendency:set_basemodels_progvarset'
    INTEGER :: jv
    INTEGER :: status
    INTEGER :: reprid_mid3d, reprid_int3d

#ifndef ICON
    CALL get_representation_id(status, 'GP_3D_MID', reprid_mid3d)
    CALL channel_halt(substr, status)
    !
    CALL get_representation_id(status, 'GP_3D_INT', reprid_int3d)
    CALL channel_halt(substr, status)
#else
    CALL get_representation_id(status, 'UCELL_3D_MID', reprid_mid3d)
    CALL channel_halt(substr, status)
    !
    CALL get_representation_id(status, 'UCELL_3D_INT', reprid_int3d)
    CALL channel_halt(substr, status)
#endif

    ! SET BASEMODEL SPECIFIC SET OF PROGNOSTIC VARIABLES
    ! EXCLUDING TRACERS
#ifdef ECHAM5
    progvar(mtend_id_t)%lexists  =.TRUE.
    progvar(mtend_id_t)%var      = t_chaobj_cpl('g1a','tm1')
    progvar(mtend_id_t)%varte    = t_chaobj_cpl('scnbuf','tte')
    progvar(mtend_id_t)%reprid   = reprid_mid3d

    progvar(mtend_id_q)%lexists  =.TRUE.
    progvar(mtend_id_q)%var      = t_chaobj_cpl('g1a','qm1')
    progvar(mtend_id_q)%varte    = t_chaobj_cpl('scnbuf','qte')
    progvar(mtend_id_q)%reprid   = reprid_mid3d

    progvar(mtend_id_xl)%lexists =.TRUE.
    progvar(mtend_id_xl)%var     = t_chaobj_cpl('g1a','xlm1')
    progvar(mtend_id_xl)%varte   = t_chaobj_cpl('scnbuf','xlte')
    progvar(mtend_id_xl)%reprid   = reprid_mid3d

    progvar(mtend_id_xi)%lexists =.TRUE.
    progvar(mtend_id_xi)%var     = t_chaobj_cpl('g1a','xim1')
    progvar(mtend_id_xi)%varte   = t_chaobj_cpl('scnbuf','xite')
    progvar(mtend_id_xi)%reprid   = reprid_mid3d

    progvar(mtend_id_u)%lexists  =.TRUE.
    progvar(mtend_id_u)%var      = t_chaobj_cpl('g2a','um1')
    progvar(mtend_id_u)%varte    = t_chaobj_cpl('scnbuf','vom')
    progvar(mtend_id_u)%reprid   = reprid_mid3d

    progvar(mtend_id_v)%lexists  =.TRUE.
    progvar(mtend_id_v)%var      = t_chaobj_cpl('g2a','vm1')
    progvar(mtend_id_v)%varte    = t_chaobj_cpl('scnbuf','vol')
    progvar(mtend_id_v)%reprid   = reprid_mid3d
#endif

#ifdef COSMO
    ! In contrast to the other basemodels, TENDENCY is just budgeting the
    ! tendencies. The original COSMO way of directly calculating the "nnew"
    ! value of a prognostic variable is not changed. Thus
    ! - new start value is always calculated by nnew + MESSy tendency
    !    (thus the progvar()%varm1ptr is set to the nnew slide AND
    !    progvar()%varteptr does only contain the additional MESSy tendencies
    ! - no closure is possible for the prognostic variables.
    ! - this also holds for the tracer tendencies
    ! To enable the simple budgeting, add_l/g routines have been equipped
    ! with a logical switch to prevent no adding of the (MESSy) tendency
    ! (progvar()%varteptr)

    progvar(mtend_id_t)%lexists =.TRUE.
    progvar(mtend_id_t)%var     = t_chaobj_cpl('COSMO','tm1')
    progvar(mtend_id_t)%varte   = t_chaobj_cpl(modstr,'TTEND')
    progvar(mtend_id_t)%reprid  = reprid_mid3d
    CALL new_channel_object(status, modstr, 'TTEND') !&
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr, 'TTEND',&
         'long_name', c='temperature tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr, 'TTEND',&
         'units',c=pvarunit(mtend_id_t))
    CALL channel_halt(substr, status)


    progvar(mtend_id_u)%lexists =.TRUE.
    progvar(mtend_id_u)%var     = t_chaobj_cpl('COSMO','um1_3d')
    progvar(mtend_id_u)%varte   = t_chaobj_cpl(modstr,'UTEND')
    progvar(mtend_id_u)%reprid  = reprid_mid3d
    CALL new_channel_object(status, modstr, 'UTEND')! &
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr, 'UTEND',&
         'long_name', c='u wind tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr, 'UTEND',&
         'units',c=pvarunit(mtend_id_u))
    CALL channel_halt(substr, status)


    progvar(mtend_id_v)%lexists =.TRUE.
    progvar(mtend_id_v)%var     = t_chaobj_cpl('COSMO','vm1_3d')
    progvar(mtend_id_v)%varte   = t_chaobj_cpl(modstr,'VTEND')
    progvar(mtend_id_v)%reprid  = reprid_mid3d

    CALL new_channel_object(status, modstr, 'VTEND')! &
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr, 'VTEND',&
         'long_name', c='v wind tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr, 'VTEND',&
         'units',c=pvarunit(mtend_id_v))
    CALL channel_halt(substr, status)

    progvar(mtend_id_w)%lexists =.TRUE.
    progvar(mtend_id_w)%var     = t_chaobj_cpl('COSMO','wm1_3d')
    progvar(mtend_id_w)%varte   = t_chaobj_cpl(modstr,'WTEND')
    progvar(mtend_id_w)%reprid  = reprid_int3d

    CALL new_channel_object(status, modstr, 'WTEND' &
         , reprid= progvar(mtend_id_w)%reprid )
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr, 'WTEND',&
         'long_name', c='w wind tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr, 'WTEND',&
         'units',c=pvarunit(mtend_id_w))
    CALL channel_halt(substr, status)

    progvar(mtend_id_pp)%lexists =.TRUE.
    progvar(mtend_id_pp)%var     = t_chaobj_cpl('COSMO','ppm1')
    progvar(mtend_id_pp)%varte   = t_chaobj_cpl(modstr,'PPTEND') ! ???
    progvar(mtend_id_pp)%reprid  = reprid_mid3d

    CALL new_channel_object(status, modstr, 'PPTEND')
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr, 'PPTEND',&
         'long_name', c='pressure perturbation tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr, 'PPTEND',&
         'units',c=pvarunit(mtend_id_pp))
    CALL channel_halt(substr, status)
#endif


#ifdef ICON
    progvar(mtend_id_t)%lexists =.TRUE.
    progvar(mtend_id_t)%var     = t_chaobj_cpl('nh_state_diag','temp')
    progvar(mtend_id_t)%varte   = t_chaobj_cpl('nh_state_diag','ddt_temp_dyn') ! ???

    progvar(mtend_id_q)%lexists =.TRUE.
    progvar(mtend_id_q)%var     = t_chaobj_cpl('tracer_gp','qv')
    progvar(mtend_id_q)%varte   = t_chaobj_cpl('tracer_gp_te','qv') ! ???
#endif
    ! GET RESPECTIVE POINTERS

    CALL main_tendency_set_domain

  END SUBROUTINE set_basemodels_progvarset

#else

  IMPLICIT NONE

#endif
end module messy_main_tendency_bi
