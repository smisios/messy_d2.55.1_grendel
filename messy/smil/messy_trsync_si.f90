! **********************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL TRSYNC 
!
! Author : Franziska Frank, DLR-IPA, 2016
!
! References: see messy_trsync.f90
!
! **********************************************************************
#include "messy_main_ppd_bi.inc"

! **********************************************************************
MODULE messy_trsync_si
! **********************************************************************
#ifdef ECHAM5

  ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      warning_bi, error_bi
  
#ifdef MESSYTENDENCY
 !tendency budget
 USE messy_main_tendency_bi,    ONLY: mtend_get_start_l,      &
                                      mtend_get_handle,       &
                                      mtend_add_l,            &
                                      mtend_register,         &
                                      mtend_id_tracer,        &
                                      mtend_id_q
#endif

  ! SMIL
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM

  ! SMCL
  USE messy_trsync

  ! GLOBAL PARAMETERS  
#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  ! MAX. NUMBER OF SYNCHRONIZED TRACER
  INTEGER, PARAMETER             :: NMAXTRSYNC = 3

  ! Array of used tracer ids
  TYPE TYPETRSYNC_ID
     INTEGER :: tr_a
     INTEGER :: tr_b
  END TYPE TYPETRSYNC_ID
  TYPE(TYPETRSYNC_ID), DIMENSION(:), POINTER :: idt_gp_trsync => NULL()

  ! Array of molarmasses
  REAL(DP), DIMENSION(:), POINTER :: my_molarmass => NULL()

  ! USED FOR NAMELIST-INPUT (ONLY INTERNALLY)
  TYPE TYPETRSYNC_IO
     CHARACTER(LEN=STRLEN_MEDIUM) :: tr_a   = ''
     CHARACTER(LEN=STRLEN_MEDIUM) :: tr_b   = ''
     INTEGER :: way = 0
  END TYPE TYPETRSYNC_IO

  ! CPL-NAMELIST PARAMETERS
  TYPE(TYPETRSYNC_IO), DIMENSION(NMAXTRSYNC), SAVE :: TRSYNC

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################
  PUBLIC :: trsync_initialize
  PUBLIC :: trsync_init_memory   ! register tendencies
  PUBLIC :: trsync_init_coupling ! set pointers for coupling to BM or other SMs
  PUBLIC :: trsync_init_tracer   ! initialize tracers
  PUBLIC :: trsync_physc         ! entry point in time loop (current vector)
  PUBLIC :: trsync_free_memory

  !--------------------------

  ! PRIVATE SUBROTINES
  !PRIVATE :: trsync_read_nml_cpl

CONTAINS

  ! ------------------------------------------------------------------------
  ! ====================================================================
  SUBROUTINE  trsync_initialize

    ! BMIL/MESSy
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,        ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'trsync_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    CALL start_message_bi(modstr, 'GLOBAL SETUP',substr)

#ifdef MESSYTENDENCY
    ! get handle for tendency treatment
    my_handle = mtend_get_handle(modstr)
#endif
    ! READ CTRL-namelist
    
    ! READ CPL-namelist
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL trsync_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    CALL p_bcast(TRSYNC%tr_a, p_io)
    CALL p_bcast(TRSYNC%tr_b, p_io)
    CALL p_bcast(TRSYNC%way, p_io)

  END SUBROUTINE trsync_initialize
  ! ====================================================================
  
  ! ====================================================================
  SUBROUTINE trsync_init_memory

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'trsync_init_memory'

#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, mtend_id_tracer)
    CALL mtend_register(my_handle, mtend_id_q)
    
#endif

  END SUBROUTINE trsync_init_memory
  ! ====================================================================
  
  ! ====================================================================
  SUBROUTINE trsync_init_coupling

    ! ------------------------------------------------------------------
    ! This soubroutine is used to set pointers
    ! (channel objects and/or tracers) for coupling to the 
    ! basemodel and to other submodels.
    ! and to transfer process-tendency-pointers from the
    ! TENDENCY-submodel to here via the mtend_request routine
    ! Specific:
    ! Couples to the submodel H2OISO as well as MECCA and CH4 
    ! respectively. It catches the HDO tendency tracer of both 
    ! respective submodels.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    ! MESSy
    USE messy_main_tracer,        ONLY: get_tracer, R_molarmass

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'trsync_init_coupling'
    INTEGER                     :: status

    INTEGER                     :: idx
    REAL(DP)                    :: molarmass_a 
    REAL(DP)                    :: molarmass_b
    INTEGER                     :: i

    INTEGER                     :: iou   ! I/O unit

    CALL start_message_bi(modstr,'TRSYNC_INIT_COUPLING',substr)  ! log-output

    ALLOCATE(idt_gp_trsync(NMAXTRSYNC))
    idt_gp_trsync(:)%tr_a = 0
    idt_gp_trsync(:)%tr_b = 0
    ALLOCATE(my_molarmass(NMAXTRSYNC))
    my_molarmass(:) = 0._dp
    molarmass_a = 0._dp
    molarmass_b = 0._dp

    DO i=1, NMAXTRSYNC

       IF (TRIM(TRSYNC(i)%tr_a) == '') CYCLE

       ! -------------

       ! Get the tracer id from tr_a
       CALL get_tracer(status, GPTRSTR, TRSYNC(i)%tr_a, idx=idx)
       CALL tracer_halt(substr, status)   ! terminate if error
       idt_gp_trsync(i)%tr_a = idx        ! save tracer id

       ! Get the tracer id from tr_b
       CALL get_tracer(status, GPTRSTR, TRSYNC(i)%tr_b, idx=idx)
       CALL tracer_halt(substr, status)   ! terminate if error
       idt_gp_trsync(i)%tr_b = idx        ! save tracer id

       ! -------------

       ! Get the molarmass from tr_a
       CALL get_tracer(status, GPTRSTR, idt_gp_trsync(i)%tr_a, &
            R_molarmass, r=molarmass_a)
       CALL tracer_halt(substr, status)   ! terminate if error

       ! Get the molarmass from tr_b
       CALL get_tracer(status, GPTRSTR, idt_gp_trsync(i)%tr_b, &
            R_molarmass, r=molarmass_b)
       CALL tracer_halt(substr, status)   ! terminate if error

       IF (.NOT. molarmass_a == molarmass_b ) THEN
          CALL error_bi('TRACER TO BE SYNCHRONIZED DO NOT HAVE'//&
               &' THE SAME MOLARMASSES!', substr)
       ELSE
          my_molarmass(i) = molarmass_a
       END IF

       ! -------------

    END DO

    CALL end_message_bi(modstr,'TRSYNC_INIT_COUPLING',substr)

  END SUBROUTINE trsync_init_coupling
  ! ====================================================================

 ! ====================================================================
  SUBROUTINE trsync_init_tracer

    ! ------------------------------------------------------------------
    ! This subroutine is used to initialise the additional
    ! water tracers according to the respective phase
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_tracer,          ONLY: get_tracer, tracer_iniflag
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev, ngpblks
    USE messy_main_data_bi,         ONLY:q, eps
    ! MESSy/SMCL
    USE messy_main_timer,         ONLY: lresume, time_step_len
!#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi, ONLY: xt, xtm1
!#endif

    USE messy_main_constants_mem, ONLY: M_air

    IMPLICIT NONE

    ! LOCAL
    INTEGER                     :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'trsync_init_tracer'
    INTEGER                     :: idt
    LOGICAL                     :: linita 
    LOGICAL                     :: linitb
    INTEGER                     :: i

    REAL(DP), DIMENSION(:,:,:), POINTER :: tr_a     => NULL() 
    REAL(DP), DIMENSION(:,:,:), POINTER :: tr_b     => NULL() 
    REAL(DP), DIMENSION(:,:,:), POINTER :: spechum  => NULL()

    REAL(DP), DIMENSION(:,:,:), POINTER :: tr_b_extra  => NULL() 

    CALL start_message_bi(modstr, 'TRSYNC tracer initialization', substr)

    ! INIT
    linita = .FALSE.
    linitb = .FALSE.

    IF (lresume) RETURN

! op_pj_20160804+ moved to here to avoid memory leak
    ! LOCAL SPACE
    ALLOCATE(tr_a(_RI_XYZ__(nproma,ngpblks,nlev)))
    ALLOCATE(tr_b(_RI_XYZ__(nproma,ngpblks,nlev)))
    ALLOCATE(spechum(_RI_XYZ__(nproma,ngpblks,nlev)))
    ALLOCATE(tr_b_extra(_RI_XYZ__(nproma,ngpblks,nlev)))
! op_pj_20160804-

    DO i=1, NMAXTRSYNC

       IF (TRIM(TRSYNC(i)%tr_a) == '') CYCLE

       ! SET POINTER TO TRACER
       CALL tracer_iniflag(status, GPTRSTR, idt_gp_trsync(i)%tr_a, lget=linita)
       CALL tracer_halt(substr, status)

       CALL tracer_iniflag(status, GPTRSTR, idt_gp_trsync(i)%tr_b, lget=linitb)
       CALL tracer_halt(substr, status)
       
       IF ((.NOT. linita).AND.(.NOT. linitb )) THEN
          CALL error_bi('NONE TRSYNC TRACER WAS INITIALISED!', substr)
       END IF

       IF (linita.AND.linitb) THEN
          !! If initializations of both tracer were carried out, only tracer
          !! with the units mol/mol_dryair will be used. Hence, the initialization
          !! of the other tracer will be overwritten.
          CALL warning_bi('BOTH TRACERS WERE INITIALIZED,'//&
            &' HOWEVER ONLY ONE CAN BE USED!'//&
!            &' => ONLY TRACER IN MOL/MOL WILL'//&
!            &' BE USED FOR INITIALIZATION!', substr)
!          linitb = .FALSE.
            &' => ONLY TRACER IN KG/KG WILL'//&
            &' BE USED FOR INITIALIZATION!', substr)
          linita = .FALSE.
       END IF

       ! ------------------------------------------------------------------
       ! Initialization via second tracer
       ! ------------------------------------------------------------------

       ! INITIALIZE
       tr_a(:,:,:) = 0.0_dp
       tr_b(:,:,:) = 0.0_dp
       spechum(:,:,:) = 0.0_dp
       tr_b_extra(:,:,:) = 0.0_dp

       ! 1. Get tracer
! Do not use start values...?? QQQQ
!#ifndef MESSYTENDENCY
!        spechum(:,:,:)  = qm1_3d(:,:,:) &
!             + qte_3d(:,:,:) * time_step_len
       idt = idt_gp_trsync(i)%tr_a
       tr_a(:,:,:) = xt(_RI_XYZN_(:,:,:,idt))
       idt = idt_gp_trsync(i)%tr_b
       tr_b(:,:,:) = xt(_RI_XYZN_(:,:,:,idt))
       spechum(:,:,:)  = q(:,:,:)   ! q == h2oisohhovap
! #else
!        CALL mtend_get_start_g(mtend_id_tracer, v0=tr_a, &
!             idt=idt_gp_trsync(i)%tr_a)
!        CALL mtend_get_start_g(mtend_id_tracer, v0=tr_b, &
!             idt=idt_gp_trsync(i)%tr_b)
!        CALL mtend_get_start_g(mtend_id_q, v0=spechum)
! #endif

       IF (.NOT. linita ) THEN
          ! init tracer a by tracer b
          ! => tr_a := convert_to_molmol(tr_b)
          
          ! 1. Convert to mol/mol
          CALL convert_unit(tr_b(:,:,:) &
               , 1                      &
               , 1                      &
               , my_molarmass(i)        &
               , spechum(:,:,:)         &
               )  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! THIS IS WORKAROUND FOR I1H2O, WHICH DOES NOT HAVE AN EXACT COUNTERPART 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF (TRSYNC(i)%tr_a == 'I1H2O') THEN
             IF (TRSYNC(1)%tr_b == 'H2OISOHDOvap') THEN
                tr_b(:,:,:) = q(:,:,:) ! correction if H2OISIHHOvap != q
                CALL convert_unit(tr_b(:,:,:) &
                     , 1                      &
                     , 1                      &
                     , my_molarmass(i)        &
                     , spechum(:,:,:)         &
                     ) 
                idt = idt_gp_trsync(1)%tr_b
                tr_b_extra(:,:,:) = xt(_RI_XYZN_(:,:,:,idt))
                CALL convert_unit(tr_b_extra(:,:,:) &
                     , 1                      &
                     , 1                      &
                     , my_molarmass(1)        &
                     , spechum(:,:,:)         &
                     ) 
                tr_b(:,:,:) = tr_b(:,:,:) - 2.0_dp*tr_b_extra(:,:,:)
             END IF
          END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! 2. Overwrite tracer
          idt = idt_gp_trsync(i)%tr_a
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! THIS IS A CRUEL WORKAROUND 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! In H2OISO istopologue tracers are declared with their
!!! factor of the atomar ratio:
!!! ratio of H atoms : 
!!! HDO / (2 * H2O + HDO) m.o.l. HDO / (2*H2O) = (0.5 * HDO) / H2O
!!!  => H2OISOHDOvap = 0.5 * HDO
!!! and          HDO = 2.0 * H2OISOHDOvap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF ( TRSYNC(i)%tr_b == 'H2OISOHDOvap' ) THEN
             xt(_RI_XYZN_(:,:,:,idt)) = 2.0_dp * tr_b(:,:,:)
          ELSE
             xt(_RI_XYZN_(:,:,:,idt)) = tr_b(:,:,:)
          END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! timestep - 1
          xtm1(_RI_XYZN_(:,:,:,idt)) = (1._DP - eps) * xt(_RI_XYZN_(:,:,:,idt))
          CALL tracer_iniflag(status, GPTRSTR, idt, lset=.TRUE.)
       END IF

       IF (.NOT. linitb) THEN
          ! init tracer b by tracer a
          ! => tr_b := convert_to_kgkg(tr_a)

          ! 1. Convert to kg/kg
          CALL convert_unit(tr_a(:,:,:) &
               , 2                      &
               , 1                      &
               , my_molarmass(i)        &
               , spechum(:,:,:)         &
               )  

          ! 2. Overwrite tracer
          idt = idt_gp_trsync(i)%tr_b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! THIS IS A CRUEL WORKAROUND 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! In H2OISO istopologue tracers are declared with their
!!! factor of the atomar ratio:
!!! ratio of H atoms : 
!!! HDO / (2 * H2O + HDO) m.o.l. HDO / (2*H2O) = (0.5 * HDO) / H2O
!!!  => H2OISOHDOvap = 0.5 * HDO
!!! and          HDO = 2.0 * H2OISOHDOvap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF ( TRSYNC(i)%tr_b == 'H2OISOHDOvap' ) THEN
             xt(_RI_XYZN_(:,:,:,idt)) = 0.5_dp * tr_a(:,:,:)
          ELSE
             xt(_RI_XYZN_(:,:,:,idt)) = tr_a(:,:,:)
          END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! timestep - 1
          xtm1(_RI_XYZN_(:,:,:,idt)) = (1._DP - eps) * xt(_RI_XYZN_(:,:,:,idt))
          CALL tracer_iniflag(status, GPTRSTR, idt, lset=.TRUE.)
       END IF
       
    END DO

    DEALLOCATE(tr_a)  ; NULLIFY(tr_a)
    DEALLOCATE(tr_b)  ; NULLIFY(tr_b)
    DEALLOCATE(spechum)  ; NULLIFY(spechum)
    DEALLOCATE(tr_b_extra)  ; NULLIFY(tr_b_extra)

    CALL end_message_bi(modstr, 'TRSYNC tracer initialization', substr)

  END SUBROUTINE trsync_init_tracer  
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE trsync_physc(flag)
    
    ! ------------------------------------------------------------------
    ! This subroutine is called within the time loop.
    ! It carries out the synchronization of the HDO tendencies from 
    ! H2OISO as well as MECCA and CH4 respectively.
    ! With flag = 1: synchronize tr_a with tr_b (=> tr_a will be overwritten)
    ! With flag = 2: synchronize tr_b with tr_a (=> tr_b will be overwritten)
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_timer,           ONLY: time_step_len
    USE messy_main_grid_def_mem_bi, ONLY: kproma, nlev, jrow
    USE messy_main_data_bi,         ONLY: qm1_3d, qte_3d

!#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, pxtm1 => qxtm1
!#endif

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*),       PARAMETER :: substr = 'trsync_physc'
    INTEGER                           :: status
    INTEGER                           :: idt
    INTEGER                           :: i

    REAL(DP), DIMENSION(:,:), POINTER :: tr_a_te     => NULL() 
    REAL(DP), DIMENSION(:,:), POINTER :: tr_b_te     => NULL() 
    REAL(DP), DIMENSION(:,:), POINTER :: tr_a        => NULL() 
    REAL(DP), DIMENSION(:,:), POINTER :: tr_b        => NULL() 
    REAL(DP), DIMENSION(:,:), POINTER :: spechum     => NULL()

    REAL(DP), DIMENSION(:,:), POINTER :: tr_b_extra  => NULL() 

!    CALL start_message_bi(modstr,'TRSYNC_PHYSC',substr)

    ! LOCAL SPACE
    ALLOCATE(tr_a_te(kproma, nlev))
    ALLOCATE(tr_b_te(kproma, nlev))
    ALLOCATE(tr_a(kproma, nlev))
    ALLOCATE(tr_b(kproma, nlev))
    ALLOCATE(tr_b_extra(kproma, nlev))
    ALLOCATE(spechum(kproma, nlev))

    tr_a_te(:,:) = 0.0_dp
    tr_b_te(:,:) = 0.0_dp
    tr_a(:,:) = 0.0_dp
    tr_b(:,:) = 0.0_dp
    spechum(:,:) = 0.0_dp

    tr_b_extra(:,:) = 0.0_dp

    DO i=1, NMAXTRSYNC

       IF (TRIM(TRSYNC(i)%tr_a) == '') CYCLE

       IF (TRSYNC(i)%way /= flag) THEN
          IF (TRSYNC(i)%way /= 0) CYCLE
       ENDIF
       
#ifndef MESSYTENDENCY
       spechum(:,:)  = qm1_3d(_RI_XYZ__(1:kproma,jrow,:)) &
            + qte_3d(_RI_XYZ__(1:kproma,jrow,:)) * time_step_len
#else
       CALL mtend_get_start_l(mtend_id_q, v0=spechum)
#endif

       idt = idt_gp_trsync(i)%tr_a
       tr_a_te(:,:) = pxtte(_RI_X_ZN_(1:kproma,:,idt))
       idt = idt_gp_trsync(i)%tr_b
       tr_b_te(:,:) = pxtte(_RI_X_ZN_(1:kproma,:,idt))

       SELECT CASE (flag)
       CASE(1)
          
          idt = idt_gp_trsync(i)%tr_a
          tr_a(:,:) = pxtm1(_RI_X_ZN_(1:kproma,:,idt)) 
!&        No tendency on tr_a, for requirements of used equation
!               + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len

#ifndef MESSYTENDENCY
          idt = idt_gp_trsync(i)%tr_b
          tr_b(:,:) = pxtm1(_RI_X_ZN_(1:kproma,:,idt))  &
               + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
#else
!!$          CALL mtend_get_start_l(mtend_id_tracer, v0=tr_b, &
!!$               idt=idt_gp_trsync(i)%tr_b)
          CALL mtend_get_start_l(idt_gp_trsync(i)%tr_b, v0=tr_b)
#endif

          ! ------------------------------------------------------------------
          ! Synchronization tr_a with tr_b (kg/kg -> mol/mol)
          ! ------------------------------------------------------------------
          ! 1. Convert tracer tendencies of tr_b to tendencies of units of tr_a_te
          ! 2. Calculate tendency difference of converted tr_b_te to tr_a_te
          !    and add it to tracer tendency

          ! 1. Convert tracer tendencies of tr_b to tendencies of units of tr_a_te
          CALL convert_unit(tr_b(1:kproma,:) &
               , 1                       &
               , 1                       &
               , my_molarmass(i)         &
               , spechum(1:kproma,:)     &
               )   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! THIS IS WORKAROUND FOR I1H2O, WHICH DOES NOT HAVE AN EXACT COUNTERPART 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF (TRSYNC(i)%tr_a == 'I1H2O') THEN
             IF (TRSYNC(1)%tr_b == 'H2OISOHDOvap') THEN
                tr_b(:,:) = spechum(1:kproma,:) ! correction if H2OISIHHOvap != q
                CALL convert_unit(tr_b(1:kproma,:) &
                     , 1                       &
                     , 1                       &
                     , my_molarmass(i)         &
                     , spechum(1:kproma,:)     &
                     ) 

#ifndef MESSYTENDENCY
! op_pj_20170904: qqq does not compile!
#else
!!$                CALL mtend_get_start_l(mtend_id_tracer, v0=tr_b_extra, &
!!$                     idt=idt_gp_trsync(1)%tr_b)  ! H2OISOHDOvap             
                CALL mtend_get_start_l(idt_gp_trsync(1)%tr_b, v0=tr_b_extra)  ! H2OISOHDOvap             
#endif

                CALL convert_unit(tr_b_extra(1:kproma,:) &
                     , 1                       &
                     , 1                       &
                     , my_molarmass(1)         &
                     , spechum(1:kproma,:)     &
                     ) 
                tr_b(:,:) = tr_b(1:kproma,:) - 2.0_dp*tr_b_extra(1:kproma,:)
             END IF
          END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! 2. Overwrite converted tracer tendency of tr_b onto tr_a_te
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! THIS IS A CRUEL WORKAROUND 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! In H2OISO istopologue tracers are declared with their
!!! factor of the atomar ratio:
!!! ratio of H atoms : 
!!! HDO / (2 * H2O + HDO) m.o.l. HDO / (2*H2O) = (0.5 * HDO) / H2O
!!!  => H2OISOHDOvap = 0.5 * HDO
!!! and          HDO = 2.0 * H2OISOHDOvap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef MESSYTENDENCY 
          idt = idt_gp_trsync(i)%tr_a
          ! tr_b_te is now in units mol/mol similar to tr_a
          IF ( TRSYNC(i)%tr_b == 'H2OISOHDOvap' ) THEN
             pxtte(_RI_X_ZN_(1:kproma,:,idt)) = ( 2.0_dp * tr_b(1:kproma,:) - tr_a(1:kproma,:) ) / time_step_len
          ELSE
             pxtte(_RI_X_ZN_(1:kproma,:,idt)) = ( tr_b(1:kproma,:) - tr_a(1:kproma,:) ) / time_step_len
          END IF
#else
          !    new tendency = converted tendency - old tendency 
          ! => sum tendency = old tendency + new tendency 
          !                 = old tendency + convert tend - old tend
          !                 = convert tendency              
          IF ( TRSYNC(i)%tr_b == 'H2OISOHDOvap' ) THEN
             tr_b_te(1:kproma,:) = ( 2.0_dp * tr_b(1:kproma,:) - tr_a(1:kproma,:) ) / time_step_len - tr_a_te(1:kproma,:)
          ELSE
             tr_b_te(1:kproma,:) = ( tr_b(1:kproma,:) - tr_a(1:kproma,:) ) / time_step_len - tr_a_te(1:kproma,:)
          END IF
          ! Add tendency of H2OISOHDOvap
#ifndef MESSYTENDENCY
! op_pj_20170904: qqq see above
#else
!!$          CALL mtend_add_l(my_handle, mtend_id_tracer, px=tr_b_te, &
!!$               idt=idt_gp_trsync(i)%tr_a)
          CALL mtend_add_l(my_handle, idt_gp_trsync(i)%tr_a, px=tr_b_te)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
          
       CASE(2)
          
          idt = idt_gp_trsync(i)%tr_b
          tr_b(:,:) = pxtm1(_RI_X_ZN_(1:kproma,:,idt)) 
!&        No tendency on tr_b, for requirements of used equation
!               + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len

#ifndef MESSYTENDENCY
          idt = idt_gp_trsync(i)%tr_a
          tr_a(:,:) = pxtm1(_RI_X_ZN_(1:kproma,:,idt))  &
               + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
#else
!!$          CALL mtend_get_start_l(mtend_id_tracer, v0=tr_a, &
!!$               idt=idt_gp_trsync(i)%tr_a)
          CALL mtend_get_start_l(idt_gp_trsync(i)%tr_a, v0=tr_a)
#endif

          ! ------------------------------------------------------------------
          ! Synchronization tr_b with tr_a (mol/mol -> kg/kg)
          ! ------------------------------------------------------------------ 
          ! 1. Convert tracer tendencies of tr_a to tendencies of units of tr_b_te
          ! 2. Calculate tendency difference of converted tr_a_te to tr_b_te
          !    and add it to tracer tendency

          ! 1. Convert tracer tendencies of tr_a to tendencies of units of tr_b_te
          CALL convert_unit(tr_a(1:kproma,:) &
               , 2                       &
               , 1                       &
               , my_molarmass(i)         &
               , spechum(1:kproma,:)     &
               )   
          ! 2. Overwrite converted tracer tendency of tr_a onto tr_b_te          
          !    and add it to tracer tendency
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! THIS IS A CRUEL WORKAROUND 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! In H2OISO istopologue tracers are declared with their
!!! factor of the atomar ratio:
!!! ratio of H atoms : 
!!! HDO / (2 * H2O + HDO) m.o.l. HDO / (2*H2O) = (0.5 * HDO) / H2O
!!!  => H2OISOHDOvap = 0.5 * HDO
!!! and          HDO = 2.0 * H2OISOHDOvap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef MESSYTENDENCY 
          idt=idt_gp_trsync(i)%tr_b
          ! tr_a is now in units kg/kg similar to tr_b
          IF ( TRSYNC(i)%tr_b == 'H2OISOHDOvap' ) THEN
             pxtte(_RI_X_ZN_(1:kproma,:,idt)) = ( 0.5_dp * tr_a(1:kproma,:) - tr_b(1:kproma,:) ) / time_step_len 
          ELSE
             pxtte(_RI_X_ZN_(1:kproma,:,idt)) = ( tr_a(1:kproma,:) - tr_b(1:kproma,:) ) / time_step_len
          END IF
#else
          !    new tendency = converted tendency - old tendency  
          ! => sum tendency = old tendency + new tendency
          !                 = old tendency + convert tend - old tend
          !                 = convert tendency              
          IF ( TRSYNC(i)%tr_b == 'H2OISOHDOvap' ) THEN
             tr_a_te(1:kproma,:) = ( 0.5_dp * tr_a(1:kproma,:) - tr_b(1:kproma,:) ) / time_step_len - tr_b_te(1:kproma,:)
          ELSE
             tr_a_te(1:kproma,:) = ( tr_a(1:kproma,:) - tr_b(1:kproma,:) ) / time_step_len - tr_b_te(1:kproma,:)
          END IF

          ! Add tendency of HDO
!!$          CALL mtend_add_l(my_handle, mtend_id_tracer, px=tr_a_te, &
!!$               idt=idt_gp_trsync(i)%tr_b)
          CALL mtend_add_l(my_handle, idt_gp_trsync(i)%tr_b, px=tr_a_te)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif             
       END SELECT

    END DO

    ! CLEAN MEMORY
    DEALLOCATE(tr_a_te)    ; NULLIFY(tr_a_te)
    DEALLOCATE(tr_b_te)    ; NULLIFY(tr_b_te)
    DEALLOCATE(tr_a)       ; NULLIFY(tr_a)
    DEALLOCATE(tr_b)       ; NULLIFY(tr_b)
    DEALLOCATE(tr_b_extra) ; NULLIFY(tr_b_extra)
    DEALLOCATE(spechum)    ; NULLIFY(spechum)

!    CALL end_message_bi(modstr,'TRSYNC_PHYSC',substr)

  END SUBROUTINE trsync_physc
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE trsync_free_memory

    IMPLICIT NONE

  END SUBROUTINE trsync_free_memory
  ! ====================================================================

  ! ####################################################################
  ! PRIVATE SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE trsync_read_nml_cpl(status, iou)
   
    ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    
    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr='trsync_read_nml_cpl'

    NAMELIST /CPL/ TRSYNC

    ! LOCAL
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE trsync_read_nml_cpl
 ! ====================================================================
#endif
! **********************************************************************
END MODULE messy_trsync_si
! **********************************************************************
