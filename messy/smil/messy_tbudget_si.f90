#include "messy_main_ppd_bi.inc"

! **********************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL TBUDGET 
!
! Author : Phoebe Graf, DLR-IPA, April 2013
!
! References: see messy_tbudget.f90
!
! **********************************************************************

! **********************************************************************
MODULE messy_tbudget_si
! **********************************************************************

  ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      error_bi, info_bi, warning_bi

  USE messy_main_tools,         ONLY: PTR_3D_ARRAY
  USE messy_main_channel,       ONLY: t_chaobj_cpl
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
#ifdef MESSYTENDENCY
 !tendency budget
 USE messy_main_tendency_bi,    ONLY: mtend_get_handle,       &
                                      mtend_get_start_l,      &
                                      mtend_add_l,            &
                                      mtend_register,           &    
                                      mtend_id_tracer
#endif

  ! SMCL
  USE messy_tbudget

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE
  SAVE

  TYPE t_dgl_io
     ! TRACER NAME, TRACER SUBNAME
     CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(2) :: tn = ''  ! new
     REAL(DP)                                   :: mm = 0.0_dp
     CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(2) :: tt = ''  ! total
     CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(2) :: tp = ''  ! prod
     ! CHANNEL OBJECT
     TYPE(t_chaobj_cpl)                         :: ol       ! loss
  END TYPE t_dgl_io

  TYPE t_dgl
     TYPE(t_dgl_io) :: io
     INTEGER :: idt_tn = 0                       ! tracer index (new diag)
     INTEGER :: idt_tt = 0                       ! tracer index (total)
     INTEGER :: idt_tp = 0                       ! tracer index (prod)
     REAL(DP), DIMENSION(:,:,:), POINTER :: pol => NULL()  ! loss rate
     INTEGER :: ix = 0   ! index of total in list
  END TYPE t_dgl

  INTEGER, PARAMETER :: NDGLMAX = 50 ! max number of equations
  INTEGER            :: NDGL    = 0  ! actual number of equations
  TYPE(t_dgl_io), DIMENSION(NDGLMAX) :: dtrac
  TYPE(t_dgl),    DIMENSION(NDGLMAX) :: xdtrac

  ! FOR LINEARISATION / RESCALING
  INTEGER :: NTOT   ! number of different totats
  INTEGER, DIMENSION(NDGLMAX) :: tot_id = 0

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  ! PUBLIC SUBROUTINES (called from messy_main_control_e5.f90)
  ! NOTE: in case you activate further entry points, make sure to call them
  !       in messy_main_control_e5.f90
  PUBLIC :: tbudget_initialize    ! initialize submodel
  PUBLIC :: tbudget_new_tracer    ! define new tracers
  PUBLIC :: tbudget_init_memory   ! define memory ...
  PUBLIC :: tbudget_init_coupling ! set pointers for coupling to BM and other SMs
  PUBLIC :: tbudget_local_end     ! entry point in time loop (current vector)

  ! PRIVATE SUBROTINES
  !PRIVATE :: tbudget_read_nml_cpl

CONTAINS

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE tbudget_initialize

    ! ------------------------------------------------------------------
    ! This subroutine is used to
!!$    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,     ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tbudget_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: jn     ! dgl counter

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

!!$    ! READ CTRL namelist
!!$    IF (p_parallel_io) THEN                  ! read only on I/O-PE
!!$       iou = find_next_free_unit(100,200)    ! find free I/O unit
!!$       CALL tbudget_read_nml_ctrl(status, iou)  ! read CTRL-namelist
!!$       ! terminate if error
!!$       IF (status /= 0) CALL error_bi('Error in reading CTRL namelist',substr)
!!$    END IF

    ! READ CPL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find next free I/O unit
       CALL tbudget_read_nml_cpl(status, iou)  ! read CPL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CPL namelist',substr)
    END IF
    ! BROADCAST CPL namleist entries from I/O-PE to ALL OTHER PEs
    DO jn = 1, NDGLMAX
       CALL p_bcast(dtrac(jn)%tn(1), p_io)
       CALL p_bcast(dtrac(jn)%tn(2), p_io)
       CALL p_bcast(dtrac(jn)%mm,    p_io)
       CALL p_bcast(dtrac(jn)%tp(1), p_io)
       CALL p_bcast(dtrac(jn)%tp(2), p_io)
       CALL p_bcast(dtrac(jn)%tt(1), p_io)
       CALL p_bcast(dtrac(jn)%tt(2), p_io)
       !
       CALL p_bcast(dtrac(jn)%ol%cha, p_io)
       CALL p_bcast(dtrac(jn)%ol%obj, p_io)
    END DO
    ! COUNT ACTUAL NUMBER OF DGLs AND SET WORKSPACE
    NDGL = 0
    DO jn = 1, NDGLMAX
       IF (TRIM(dtrac(jn)%tn(1)) == '') CYCLE
       IF (TRIM(dtrac(jn)%tp(1)) == '') CYCLE
       IF (TRIM(dtrac(jn)%tt(1)) == '') CYCLE
       IF (TRIM(dtrac(jn)%ol%cha) == '') CYCLE
       IF (TRIM(dtrac(jn)%ol%obj) == '') CYCLE
       NDGL = NDGL + 1
       xdtrac(NDGL)%io%tn(1)   = dtrac(jn)%tn(1)
       xdtrac(NDGL)%io%tn(2)   = dtrac(jn)%tn(2)
       xdtrac(NDGL)%io%mm      = dtrac(jn)%mm
       xdtrac(NDGL)%io%tp(1)   = dtrac(jn)%tp(1)
       xdtrac(NDGL)%io%tp(2)   = dtrac(jn)%tp(2)
       xdtrac(NDGL)%io%tt(1)   = dtrac(jn)%tt(1)
       xdtrac(NDGL)%io%tt(2)   = dtrac(jn)%tt(2)
       xdtrac(NDGL)%io%ol%cha  = dtrac(jn)%ol%cha
       xdtrac(NDGL)%io% ol%obj = dtrac(jn)%ol%obj
    END DO

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif

    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

  END SUBROUTINE tbudget_initialize
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE tbudget_new_tracer

    ! ------------------------------------------------------------------
    ! This subroutine is used to define new tracers. See
    ! http://www.atmos-chem-phys.net/8/1677   (including supplement !)
    ! for full documentation.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    ! MESSy
    USE messy_main_tracer,        ONLY: new_tracer, set_tracer &
                                      , R_molarmass ! ,ON, OFF
    USE messy_main_constants_mem, ONLY: MO

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tbudget_new_tracer'
    INTEGER                     :: status
    INTEGER                     :: jn     ! dgl counter

    CALL start_message_bi(modstr,'TRACER DEFINITION',substr)  ! log-output

    DO jn = 1, NDGL
       CALL new_tracer(status, GPTRSTR, TRIM(xdtrac(jn)%io%tn(1)) &
            , modstr, idx=xdtrac(jn)%idt_tn, subname=TRIM(xdtrac(jn)%io%tn(2)) &
            , unit='mol/mol' )
       CALL tracer_halt(substr, status)   ! terminate if error
       CALL set_tracer(status, GPTRSTR, xdtrac(jn)%idt_tn &
         , R_molarmass, r=xdtrac(jn)%io%mm)
       CALL tracer_halt(substr, status)   ! terminate if error
    END DO

    CALL end_message_bi(modstr,'TRACER DEFINITION',substr)  ! log-output

  END SUBROUTINE tbudget_new_tracer
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE tbudget_init_memory

    IMPLICIT NONE

#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, mtend_id_tracer)
#endif

  END SUBROUTINE tbudget_init_memory
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE tbudget_init_coupling

    ! ------------------------------------------------------------------
    ! This soubroutine is used to set pointers
    ! (channel objects and/or tracers) for coupling to the 
    ! basemodel and to other submodels.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: get_channel_object
    !
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_tracer,        ONLY: get_tracer
    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tbudget_init_coupling'
    INTEGER                     :: status
    INTEGER                     :: jn     ! dgl counter
    INTEGER                     :: i
    LOGICAL                     :: lex

    CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output

    DO jn = 1, NDGL
       ! total
       CALL get_tracer(status, GPTRSTR, TRIM(xdtrac(jn)%io%tt(1)) &
            ,TRIM(xdtrac(jn)%io%tt(2)), xdtrac(jn)%idt_tt)
       CALL tracer_halt(substr, status)   ! terminate if error

       ! prod
       CALL get_tracer(status, GPTRSTR, TRIM(xdtrac(jn)%io%tp(1)) &
            ,TRIM(xdtrac(jn)%io%tp(2)), xdtrac(jn)%idt_tp)
       CALL tracer_halt(substr, status)   ! terminate if error

       ! loss
       CALL get_channel_object(status, TRIM(xdtrac(jn)%io%ol%cha) &
            , TRIM(xdtrac(jn)%io%ol%obj), p3=xdtrac(jn)%pol)
       CALL channel_halt(substr, status)
    END DO

    ! save indices for linearisation / rescaling
    NTOT = 0
    DO jn = 1, NDGL
       lex = .FALSE.
       DO i=1, NTOT
          IF (xdtrac(jn)%idt_tt == tot_id(i)) THEN
             lex = .TRUE.  ! exists already in list
             EXIT
          END IF
       END DO
       IF (.NOT. lex) THEN
          NTOT = NTOT + 1
          tot_id(NTOT) = xdtrac(jn)%idt_tt
          xdtrac(jn)%ix = NTOT
       ELSE
          xdtrac(jn)%ix = i
       END IF
    END DO

    CALL end_message_bi(modstr,'COUPLING',substr)  ! log-output

  END SUBROUTINE tbudget_init_coupling
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE tbudget_local_end    

#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi, ONLY: qxtm1, qxtte
#else
    USE messy_main_tracer_mem_bi, ONLY: qxtte
#endif
    USE messy_main_grid_def_mem_bi, ONLY: kproma, nlev, jrow
    USE messy_main_timer,         ONLY: time_step_len

    IMPLICIT NONE
    INTRINSIC :: EPSILON

    ! LOCAL
    INTEGER :: jn ! loop counter
    INTEGER :: idt   ! tracer ID
    ! start value for new diagnostic (tagged) tracer
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ztn0
    ! start value for total budget tracer
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ztt0
    ! additonal tendency for new diagnostic (tagged) tracer
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ztnte
    ! correction factor to force sum of tagged tracers being equal to
    ! total budget tracer
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: acorr

    ! for each tagging DGL
    ALLOCATE(ztnte(kproma,nlev,NDGL))
    ztnte(:,:,:) = 0.0_dp
    ALLOCATE(ztn0(kproma,nlev,NDGL))
    ztn0(:,:,:) = 0.0_dp

    ! for each total budget tracer
    ALLOCATE(acorr(kproma,nlev,NTOT))
    acorr(:,:,:) = 0.0_dp
    ALLOCATE(ztt0(kproma,nlev,NTOT))
    ztt0(:,:,:) = 0.0_dp

    loop1: DO jn=1, NDGL

#ifndef MESSYTENDENCY
       idt = xdtrac(jn)%idt_tn
       ztn0(:,:,jn) = qxtm1(_RI_X_ZN_(1:kproma,:,idt)) + &
                      qxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len

       idt = xdtrac(jn)%idt_tt
       ztt0(:,:,xdtrac(jn)%ix) = qxtm1(_RI_X_ZN_(1:kproma,:,idt)) + &
                                 qxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len  
#else
!!$       CALL mtend_get_start_l(mtend_id_tracer &
!!$            , v0=ztn0(:,:,jn), idt=xdtrac(jn)%idt_tn)
!!$       CALL mtend_get_start_l(mtend_id_tracer &
!!$            , v0=ztt0(:,:,xdtrac(jn)%ix), idt=xdtrac(jn)%idt_tt)
       CALL mtend_get_start_l(xdtrac(jn)%idt_tn, v0=ztn0(:,:,jn))
       CALL mtend_get_start_l(xdtrac(jn)%idt_tt, v0=ztt0(:,:,xdtrac(jn)%ix))
#endif

       idt = xdtrac(jn)%idt_tp
       CALL dgl_step(ztnte(:,:,jn)             &  ! OUT: new tendency
            , qxtte(_RI_X_ZN_(1:kproma,:,idt))         &  ! mol/mol/s (prod)
            , xdtrac(jn)%pol(_RI_XYZ__(1:kproma,jrow,:)) &  ! mol/mol/s (loss) ! < 0 !!!
            , ztn0(:,:,jn)                     &  ! mol/mol   (tagged tracer)
            , ztt0(:,:,xdtrac(jn)%ix)          &  ! mol/mol   (total budget)
            )

       ! summation of different components for each total
       ! %ix is number in list
       acorr(:,:,xdtrac(jn)%ix) = acorr(:,:,xdtrac(jn)%ix) + ztn0(:,:,jn) &
            + ztnte(:,:,jn) * time_step_len

    END DO loop1

    ! calculate correction factor as quotient of total tracer (at t+1) to
    ! sum of tagged tracers (at t+1)
    WHERE (acorr(:,:,:) >= EPSILON(0.0_dp))
       acorr(:,:,:) = ztt0(:,:,:)/acorr(:,:,:)
    ELSEWHERE
       acorr(:,:,:) = 1.0_dp
    END WHERE

    ! s + t*dt =!= (f1 + t1*dt) + (f2 + t2*dt) + ...
    ! acorr = a := (s+t*dt) / ( (f1+t1*dt) + (f2+t2*dt) + ... )
    ! (f1+t1*dt)' := a*(f1+t1*dt)
    ! (f2+t2*dt)' := a*(f2+t2*dt)
    ! ...
    ! => (f1+t1*dt)' + (f2+t2*dt)' + ... = (s+t*dt)   O.K.
    ! do not change tracers, rather adjust tendencies
    ! cond.: (f1' =!= f1) AND (f2' =!= f2) AND ...
    ! => f1 + t1'*dt = a*f1 + a*t1*dt ; f2 + t2'*dt = a*f2 + a*t2*dt ; ...
    ! => t1' = f1*(a-1)/dt + a*t1  ; t2' = f2*(a-1)/dt + a*t2
    loop2: DO jn = 1, NDGL
       ztnte(:,:,jn) = acorr(:,:,xdtrac(jn)%ix)*ztnte(:,:,jn) + &
            (acorr(:,:,xdtrac(jn)%ix) - 1.0_dp)*ztn0(:,:,jn)/time_step_len
    END DO loop2

    loop3: DO jn = 1, NDGL

       idt = xdtrac(jn)%idt_tn
#ifndef MESSYTENDENCY
       qxtte(_RI_X_ZN_(1:kproma,:,idt)) = qxtte(_RI_X_ZN_(1:kproma,:,idt)) + ztnte(:,:,jn)
#else
!!$       CALL mtend_add_l(my_handle, mtend_id_tracer &
!!$            , px= ztnte(:,:,jn), idt=idt)
       CALL mtend_add_l(my_handle, idt, px= ztnte(:,:,jn))
#endif

    END DO loop3

    DEALLOCATE(ztnte)
    DEALLOCATE(ztn0)
    DEALLOCATE(acorr)
    DEALLOCATE(ztt0)

  END SUBROUTINE tbudget_local_end
  ! ====================================================================

  ! ####################################################################
  ! PRIVATE SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE tbudget_read_nml_cpl(status, iou)
   
    ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    
    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ dtrac

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='tbudget_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! ### ADD HERE DIAGNOSTIC OUTPUT FOR LOG-FILE

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE tbudget_read_nml_cpl
  ! ====================================================================

! **********************************************************************
END MODULE messy_tbudget_si
! **********************************************************************
