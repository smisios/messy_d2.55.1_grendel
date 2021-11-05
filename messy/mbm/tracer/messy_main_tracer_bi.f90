! **********************************************************************+
MODULE messy_main_tracer_bi
! **********************************************************************+

  ! THIS MODULE PROVIDES THE BML-SPECIFIC MAIN ENTRY POINTS FOR THE
  ! GENERIC MESSy-SUBMODEL 'TRACER'

  ! BMIL
  USE messy_main_tracer_mem_bi
  USE messy_main_tools_bi,        ONLY: start_message_si, end_message_si
  USE messy_main_tracer_tools_bi, ONLY: tracer_halt
  ! SMCL
  USE messy_main_tracer

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: S1TRSTR

  ! PUBLIC SUBROUTINES CALLED FROM BM(I)L
  PUBLIC :: main_tracer_initialize    ! CALLED FROM tracer_bml
  PUBLIC :: main_tracer_new_tracer    ! CALLED FROM tracer_bml
  PUBLIC :: main_tracer_init_memory   ! CALLED FROM tracer_bml
  PUBLIC :: main_tracer_init_coupling ! CALLED FROM tracer_bml
  PUBLIC :: main_tracer_init_tracer   ! CALLED FROM tracer_bml
  PUBLIC :: main_tracer_global_start  ! CALLED FROM tracer_bml
  !
  ! SPECIAL: CALLED FROM BML, AFTER ADVECTION
  PUBLIC :: main_tracer_beforeadv     ! CALLED FROM tracer_bml
  PUBLIC :: main_tracer_afteradv      ! CALLED FROM tracer_bml
  !
  PUBLIC :: main_tracer_local_start   ! CALLED FROM tracer_bml
  PUBLIC :: main_tracer_global_end    ! CALLED FROM tracer_bml
  PUBLIC :: main_tracer_free_memory   ! CALLED FROM tracer_bml
  !
  ! CONVERSION ROUTINES TO BE USED IN SUBMODEL SMIL
  PUBLIC :: main_tracer_fconv_loc    ! CONVERT FAMILIES 1 <-> TRACERS
  PUBLIC :: main_tracer_fconv_glb    ! CONVERT FAMILIES 1 <-> TRACERS

CONTAINS

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_initialize

    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT

    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_io, p_bcast, finish
    USE messy_main_tools_bi,         ONLY: find_next_free_unit
    ! SMIL
    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_initialize
    USE messy_main_tracer_pdef_bi,   ONLY: main_tracer_pdef_initialize

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL main_tracer_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL finish(substr, &
            'main_tracer_read_nml_ctrl reported an error')
    END IF
    CALL p_bcast(l_family, p_io)
    CALL p_bcast(l_pdef, p_io)

    IF (l_family) CALL main_tracer_family_initialize

    IF (l_pdef) CALL main_tracer_pdef_initialize

  END SUBROUTINE main_tracer_initialize
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_new_tracer(flag)

    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, finish
    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_new_tracer

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_new_tracer'
    INTEGER                     :: status

    SELECT CASE(flag)
    CASE(1)
       !
       ! DEFINE NEW TRACER SET
       CALL start_message_si(modstr,'SETUP TRACER SETS',substr)
       !
       CALL new_tracer_set(status, S1TRSTR, .TRUE.)
       CALL tracer_halt(substr, status)
       !
       CALL end_message_si(modstr,'SETUP TRACER SETS',substr)
       !
    CASE (2)
       !
       IF (l_family) CALL main_tracer_family_new_tracer
       !
    CASE(3)
       !
       ! op_pj_20160823+
       CALL start_message_si(modstr,'OVERWRITE TRACER PROPERTIES',substr)
       CALL set_tracer_properties(status)
       CALl tracer_halt(substr, status)
       CALL end_message_si(modstr,'OVERWRITE TRACER PROPERTIES',substr)
       ! op_pj_20160823-
       !
       ! DIAGNOSTIC OUTPUT
       CALL start_message_si(modstr,'SHOW TRACER SETS',substr)
       !
       IF (p_parallel_io) CALL print_tracer_set
       !
       CALL end_message_si(modstr,'SHOW TRACER SETS',substr)
       !
    CASE DEFAULT
       !
       CALL finish(substr, 'UNKNOWN FLAG !')
       !
    END SELECT

  END SUBROUTINE main_tracer_new_tracer
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_init_memory(flag)

    ! BMIL
    USE messy_main_mpi_bi,           ONLY: finish
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, nlev, ngpblks
    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_init_mem
    USE messy_main_tracer_pdef_bi,   ONLY: main_tracer_pdef_init_mem

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)   :: flag
    INTEGER, DIMENSION(3) :: dims

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_init_memory'
    INTEGER :: status

    SELECT CASE(flag)
    CASE (1)
       !
       CALL start_message_si(modstr,'SETUP TRACER MEMORY',substr)
       !
       dims(1) = nproma
       dims(2) = nlev
       dims(3) = ngpblks
       !
       CALL setup_tracer_set(status, S1TRSTR, dims, 1, .FALSE., .TRUE.)
       CALL tracer_halt(substr, status)
       !
       ! SETUP POINTERS TO GRIDPOINT TRACER MEMORY AND TRACER INFORMATION
       CALL get_tracer_set(status, S1TRSTR, trlist_s1, ti_s1, ntrac_s1 &
         , xt=pxt)
       CALL tracer_halt(substr, status)
       !
       ! INITALIZE POINTER (messy_main_tracer_mem_bi.90)
       xt => pxt(:,:,:,:,1)
       !
       CALL end_message_si(modstr,'SETUP TRACER MEMORY',substr)
       !
    CASE (2)
       !
       CALL start_message_si(modstr,'TRACER MEMORY CHANNEL COUPLING',substr)
       !
       ! ###+ 
       ! ### ADD HERE: COUPLE TRACER MEMORY TO 
       ! ###           MESSy DATA TRANSFER / EXPORT INTERFACE
       ! ###-
       !
       IF (l_pdef) CALL main_tracer_pdef_init_mem
       !
       ! setting meta information of family-members to fraction
       IF (l_family) CALL main_tracer_family_init_mem
       !
       CALL end_message_si(modstr,'TRACER MEMORY CHANNEL COUPLING',substr)
       !
    CASE DEFAULT
       !
       CALL finish(substr, 'UNKNOWN FLAG !')
       !
    END SELECT

  END SUBROUTINE main_tracer_init_memory
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_init_coupling

    ! BMIL
    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_init_cpl
    USE messy_main_tracer_pdef_bi,   ONLY: main_tracer_pdef_init_cpl

    IMPLICIT NONE

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_init_coupling'

    IF (l_family) CALL main_tracer_family_init_cpl

    IF (l_pdef) CALL main_tracer_pdef_init_cpl

  END SUBROUTINE main_tracer_init_coupling
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_init_tracer(flag)

    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, finish
!!$    USE messy_main_data_bi,          ONLY: eps, lresume
    ! SMCL
    USE messy_main_tracer,           ONLY: NSETID, STRLEN_TRSET &
                                         , get_tracer_set, R_vini

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
!!$    REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: zpxtm1 => NULL()
    LOGICAL :: l_init
    INTEGER :: i, ntrac

    SELECT CASE(flag)
    CASE(1)
       !
       CALL start_message_si(modstr,'CHECK TRACER INIT FROM RERUN',substr)
       !
       WRITE(*,*) '[THIS FEATURE IS NOT IMPLEMENTED IN THIS BOX MODEL]'
       !
! op_pj_20100114+
    CASE(4)
       !
       CALL start_message_si(modstr,'INITIALISE TRACER VIA TRACER_INIT',substr)
       WRITE(*,*) '[THIS FEATURE IS NOT IMPLEMENTED IN THIS BOX MODEL]'
       !CALL tracer_init
       CALL end_message_si(modstr,'INITIALISE TRACER VIA TRACER_INIT',substr)
       RETURN
       !
! op_pj_20100114-
    CASE(2)
       !
       CALL start_message_si(modstr,'CHECK TRACER INIT BY TRACER_INIT',substr)
       !
    CASE(3)
       !
       CALL start_message_si(modstr,'DIAGNOSE TRACER INITIALIZATION',substr)
       !
    CASE DEFAULT
       !
       CALL finish(substr, 'UNKNOWN FLAG !')
       !
    END SELECT

    set_loop: DO i=1, NSETID
       !
       CALL get_tracer_set(status, i, setname=setname &
            , trlist=til, xt=zpxt &
!!$            , trlist=til, xt=zpxt, xtm1=zpxtm1 &
            , ntrac=ntrac, l_init=l_init)
       CALL tracer_halt(substr, status)

       IF (flag < 3) THEN
          IF (p_parallel_io) THEN
             WRITE(*,*) '*** TRACER SET '//TRIM(setname)//': ',&
                  ntrac,' TRACERS'
          END IF
       END IF

       IF (.NOT. l_init) CYCLE
       IF (ntrac <= 0) CYCLE

       !
       tracer_loop: DO
          IF (.NOT.ASSOCIATED(til)) EXIT

          SELECT CASE(flag)
          CASE(1)
             !
             ! ###+ 
             ! ###  ADD HERE:
             ! ###  CHECK IF TRACER HAS BEEN INITIALIZED FROM RESTART FILES:
             ! ###  MESSy DATA TRANSFER AND EXPORT INTERFACE OBJECTS
             ! ###  FOR X and X_m1 have restart_read = .true.
             ! ###  THIS IS NOT IMPLEMENTED FOR THIS BOX MODEL!
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
                        ') WILL BE SET TO ', til%info%meta%cask_r(R_vini)
                END IF
                zpxt(:,:,til%info%ident%idx,:,:) = til%info%meta%cask_r(R_vini)
!!$                IF (lresume) &
!!$                     zpxtm1(:,:,til%info%ident%idx,:,:) = (1._DP -eps) &
!!$                     * zpxt(:,:,til%info%ident%idx,:,:) 
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
       CALL end_message_si(modstr,'CHECK TRACER INIT FROM RERUN',substr)
       !
    CASE(2)
       !
       CALL end_message_si(modstr,'CHECK TRACER INIT BY TRACER_INIT',substr)
       !
    CASE(3)
       !
       ! --- DIAGNOSTIC OUTPUT
       !
       IF (p_parallel_io) CALL print_tracer_set_val
       !
       CALL end_message_si(modstr,'DIAGNOSE TRACER INITIALIZATION',substr)
       !
    CASE DEFAULT
       !
       CALL finish(substr, 'UNKNOWN FLAG !')
       !
    END SELECT

  END SUBROUTINE main_tracer_init_tracer
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_global_start(flag)

    USE messy_main_tracer_pdef_bi,   ONLY: main_tracer_pdef_global_start
! um_ak_20081030+
   ! USE messy_main_tracer_family_bi, ONLY: main_tracer_family_beforeadv
! um_ak_20081030-

    IMPLICIT NONE

    !  I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_global_start'

    SELECT CASE(flag)
       !
    CASE(1)
       !
       IF (l_pdef) CALL main_tracer_pdef_global_start
       !
    CASE(2)
       !
       ! TYPE-1: t2f
       ! TYPE-2: summation (GPTRSTR)
       !IF (l_family) CALL main_tracer_family_beforeadv ! um_ak_20081030
       !
    END SELECT

  END SUBROUTINE main_tracer_global_start
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! um_ak_20080709+
  SUBROUTINE main_tracer_beforeadv

    USE messy_main_tracer_family_bi,   ONLY: main_tracer_family_beforeadv
    
    IMPLICIT NONE

       ! TYPE-1: t2f
       ! TYPE-2: summation (GPTRSTR)
       IF (l_family) CALL main_tracer_family_beforeadv

  END SUBROUTINE main_tracer_beforeadv
  ! um_ak_20080709-
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

    ! BMIL
    USE messy_main_grid_def_mem_bi,    ONLY: jrow

    IMPLICIT NONE

    ! SET POINTERS FOR WITHIN LOCAL LOOP

    ! GRIDPOINT TRACERS
    qxt    => xt(:,:,:,jrow)
!!$    qxtte  => xtte(:,:,:,jrow)
!!$    qxtm1  => xtm1(:,:,:,jrow)
!!$    qxtf   => xtf(:,:,:,jrow)
    
  END SUBROUTINE main_tracer_local_start
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_global_end

    ! BMIL
    USE messy_main_grid_def_bi,    ONLY: grmass
    ! SMIL
    USE messy_main_tracer_pdef_bi, ONLY: main_tracer_pdef_global_end
    ! SMCL
    USE messy_main_tracer_pdef,    ONLY: tracpdef_airmass

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_global_end'

    IF (l_pdef) THEN
       CALL tracpdef_airmass(S1TRSTR, grmass)
       CALL main_tracer_pdef_global_end
    END IF

  END SUBROUTINE main_tracer_global_end
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_free_memory

    USE messy_main_tracer_pdef_bi,   ONLY: main_tracer_pdef_free_mem
    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_free_mem

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_free_memory'
    INTEGER                     :: status

    CALL start_message_si(modstr,'FREE TRACER MEMORY',substr)

    CALL clean_tracer_set(status, S1TRSTR)
    CALL tracer_halt(substr, status)

    IF (l_pdef) CALL main_tracer_pdef_free_mem

    IF (l_family) CALL main_tracer_family_free_mem

    CALL end_message_si(modstr,'FREE TRACER MEMORY',substr)

  END SUBROUTINE main_tracer_free_memory
  ! -------------------------------------------------------------------

 ! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_fconv_loc(direction, callstr, TRSETSTR)

    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_grid_def_mem_bi,  ONLY: jrow, kproma
#if defined(MBM_TRACER)
    USE messy_main_data_bi,          ONLY: time_step_len
#else
    USE messy_main_timer,            ONLY: time_step_len
#endif
    !
    ! SMCL
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

    ! BMIL
    USE messy_main_mpi_bi,          ONLY: p_pe
    USE messy_main_grid_def_mem_bi, ONLY: jrow, nproma, npromz &
                                        , ngpblks
    USE messy_main_data_bi,         ONLY: time_step_len
    !
    ! SMCL
    USE messy_main_tracer_family_bi, ONLY: tracfamily_1_f2t, tracfamily_1_t2f &
                                         , tracfamily_2_sum, tracfamily_2_rsc

    IMPLICIT NONE

    INTRINSIC :: TRIM

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
          IF (jjrow == ngpblks) THEN
             jp = npromz
          ELSE
             jp = nproma
          END IF
          CALL tracfamily_1_f2t(status, callstr, p_pe, TRSETSTR, &
               time_step_len, jjrow, jp)
          CALL tracer_halt(substr, status)
       END DO
       !
    CASE ('t2f','T2f','t2F','T2F')
       !
       DO jjrow = 1, ngpblks
          IF (jjrow == ngpblks) THEN
             jp = npromz
          ELSE
             jp = nproma
          END IF
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

! **********************************************************************+
END MODULE messy_main_tracer_bi
! **********************************************************************+
